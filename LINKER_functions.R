run_linker<-function(lognorm_est_counts, protein_filtered_idx, lincs_filtered_idx, NrModules, module_summary, 
                     Lambda=0.0001, alpha=1-1e-06,
                     pmax=10, mode="LASSO", used_method="LINKER",
                     NrCores=30, corrClustNrIter=21,
                     Nr_bootstraps=1)
{
  
  # Creating the parameters structure
  Parameters <- list(Lambda=Lambda,pmax=pmax,alpha=alpha, mode=mode, used_method=used_method)
  sample_size<-dim(lognorm_est_counts)[2]
  test_size<-round(0.1*sample_size)
  train_size<-round(0.7*sample_size)
  EvaluateTestSet<-list()
  bootstrap_modules<-list()
  bootstrap_results<-list()
  test_samples<-sample(1:sample_size, test_size, replace=F)
  train_val_samples<-setdiff(1:sample_size, test_samples)
  
  Regulator_data_test = t(scale(t(lognorm_est_counts[lincs_filtered_idx,test_samples])))
  MA_matrix_Var_test = t(scale(t(lognorm_est_counts[protein_filtered_idx,test_samples])))
  
  Regulator_data_train_val = t(scale(t(lognorm_est_counts[lincs_filtered_idx,train_val_samples])))
  MA_matrix_Var_train_val = t(scale(t(lognorm_est_counts[protein_filtered_idx,train_val_samples])))
  
  for(boost_idx in 1:Nr_bootstraps)
  {
    
    train_samples<-sample(train_val_samples, train_size, replace=F)
    validation_samples<-setdiff(train_val_samples, train_samples)
    
    Regulator_data_train = t(scale(t(lognorm_est_counts[lincs_filtered_idx,train_samples])))    
    Regulator_data_validation = t(scale(t(lognorm_est_counts[lincs_filtered_idx,validation_samples])))
    
    MA_matrix_Var_train = t(scale(t(lognorm_est_counts[protein_filtered_idx,train_samples])))
    MA_matrix_Var_validation = t(scale(t(lognorm_est_counts[protein_filtered_idx,validation_samples])))
    
    LINKERinit<-LINKER_init(MA_matrix_Var = MA_matrix_Var_train, RegulatorData = Regulator_data_train, NrModules = NrModules, NrCores=NrCores, corrClustNrIter=corrClustNrIter, Parameters = Parameters )
    
    
    tmp<-LINKER_corrClust(LINKERinit)
    bootstrap_results[[boost_idx]]<-tmp
    
    EvaluateTestSet[[boost_idx]]<-LINKER_EvaluateTestSet(bootstrap_results[[boost_idx]],MA_matrix_Var_validation,Regulator_data_validation, used_method=Parameters$used_method)
    
    printf("Bootstrap %d, NrModules %d:\n", boost_idx, bootstrap_results[[boost_idx]]$NrModules)
    
    #print(EvaluateTestSet[[boost_idx]])
    
    print(apply(EvaluateTestSet[[boost_idx]], 2, mean))
  }
  
  
  if(module_summary==0){
    return(list(bootstrapResults=bootstrap_results,bootstrapTestStats=EvaluateTestSet))
  }
  
  # bootstrap_modules[[boost_idx]]<-list()
  # for(i in 1:bootstrap_results[[boost_idx]]$NrModules){
  #   
  #   module_genes<-bootstrap_results[[boost_idx]]$AllGenes[which(bootstrap_results[[boost_idx]]$ModuleMembership==i)]
  #   module_gene_description_list<-sapply(module_genes, function(x) strsplit(x, "\\|"))
  #   bootstrap_modules[[boost_idx]][[i]]<-sapply(module_gene_description_list, function(x) x[[6]])
  # }
  
  #corrClust_results<-list(bootstrap_modules, EvaluateTestSet, bootstrap_results)
  
  bootstrap_matrix<-LINKER_filterModuleBootstraps(bootstrap_results,EvaluateTestSet)
  
  consensus_modules<-LINKER_consensusClustering_IPC(bootstrap_matrix)
  
  regulators<-LINKER_compute_lincRegulators(MA_matrix_Var_train_val, Regulator_data_train_val, consensus_modules, Parameters)
  
  testResults<-LINKER_EvaluateTestSet(regulators,MA_matrix_Var_test,Regulator_data_test, used_method=Parameters$used_method)
  
  printf("NrModules %d:\n",nrow(testResults))
  print(apply(testResults, 2, mean))
  
  return(list(consensusTestStats=testResults, consensusResults=regulators, bootstrapResults=bootstrap_results,bootstrapTestStats=EvaluateTestSet))

}

LINKER_init <- function(MA_matrix_Var, RegulatorData, NrModules, NrCores=30, corrClustNrIter=21, Parameters) {
  
  if (nrow(MA_matrix_Var)>NrModules){
    
    #Set initialization of centroids
    #rnd_cent<-1:NrModules
    
    #Random initailization of the centroids
    #rnd_cent<-runif(NrModules, min = 1, max = nrow(MA_matrix_Var))
    
    #K-means++-style initialization
    rnd_cent<-matrix()
    #First one at random
    rnd_cent<-(MA_matrix_Var[sample(1:nrow(MA_matrix_Var), 1),])
    #the second value is the least correlated from the first one
    corr_dists<-apply( MA_matrix_Var, 1, function(x) ( (x%*%rnd_cent)/(length(x)-1) )^2 )
    rnd_cent<-cbind(as.matrix(rnd_cent),as.matrix(MA_matrix_Var[which(corr_dists==min(corr_dists))[1],]))
    rnd_cent<-t(rnd_cent)
    #compute the rest of the centers
    Px<-numeric(length = nrow(MA_matrix_Var))
    for(center_idx in 3:NrModules){
      
      for(i in 1: nrow(MA_matrix_Var)){
        gene<-MA_matrix_Var[i,]
        corr_dists<-apply( rnd_cent, 1, function(x) ((x%*%gene)/(length(x)-1) )^2 )
        Px[i]<-2-max(corr_dists)
      }
      if(sum(is.finite(Px))!= length(Px)){
        print("asdf")
      }
      rnd_cent<-rbind(rnd_cent,MA_matrix_Var[sample(1:nrow(MA_matrix_Var), 1, prob = Px),])
    }
    
    
    
    
    ModuleVectors<-rnd_cent
    
    Data<-MA_matrix_Var
    Clusters<-numeric()
    for(jj in 1:5){
      for (i in 1:nrow(MA_matrix_Var)){
        CurrentGeneVector = Data[i,,drop=FALSE]
        Correlations = abs(cor(t(CurrentGeneVector),t(ModuleVectors)))
        Correlations = (cor(t(CurrentGeneVector),t(ModuleVectors)))
        corr = data.matrix(Correlations,rownames.force = NA)
        #corr[is.na(corr)] <- -100000
        MaxCorrelation = max(corr,na.rm=TRUE)
        #print(MaxCorrelation)
        #MaxPosition = which(corr == max(corr,na.rm=TRUE))
        MaxPosition = which(signif(corr,digits=7) == signif(MaxCorrelation,digits=7))
        #print(MaxPosition)
        MaxPosition = MaxPosition[1] # this is new, to avoid two different reassignements
        
        Clusters[i] = MaxPosition
      }
      #print(rle(Clusters)$length)
      ClusterIDs<-unique(Clusters)
      for (idx in 1:length(ClusterIDs)){
        genesInModule<-which(Clusters == ClusterIDs[idx])
        cx <- MA_matrix_Var[genesInModule, ]
        if(length(genesInModule) > 1){
          if(Parameters$used_method=="LINKER"){
            clusterSVD <- svd(cx)
            y<-clusterSVD$v[,1]
            if(cor(apply(cx,2,mean), y)<0){
              y<- -y
            }
            #ModuleVectors[idx,]<-scale(y, center = FALSE)
            
          }
          else{
            ModuleVectors[idx,]<-apply(cx,2,mean)
          }
        }
        else{
          ModuleVectors[idx,]<-cx
        }
        
      }
    }
    
  } else {
    stop("The number of modules is too large compared to the total number of genes.")
  }
  ModuleMembership<-as.numeric(Clusters)
  names(ModuleMembership) <- rownames(MA_matrix_Var)
  
  return(list(MA_matrix_Var=MA_matrix_Var,RegulatorData=RegulatorData,ModuleMembership=ModuleMembership,Parameters=Parameters, NrCores=NrCores, corrClustNrIter=corrClustNrIter))
  
}


LINKER_corrClust <- function(LINKERinit){
  
  
  NrIterations<-LINKERinit$corrClustNrIter
  
  
  if (nrow(LINKERinit$RegulatorData)==1){
    stop('Only one driver is detected. More than one driver is needed.\n')
  }
  #cat('Running on',length(rownames(LINKERinit$MA_matrix_Var)),'transcripts/genes and',length(colnames(LINKERinit$MA_matrix_Var)),'samples.\n')
  
  
  Data<-LINKERinit$MA_matrix_Var
  Clusters<-LINKERinit$ModuleMembership
  RegulatorData<-LINKERinit$RegulatorData
  Parameters<-LINKERinit$Parameters
  NrCores<-LINKERinit$NrCores
  
  # this will register nr of cores/threads, keep this here so the user can decide how many cores based on their hardware.
  registerDoParallel(cores=NrCores)
  ptm1 <- proc.time()
  
  RegulatorData_rownames=rownames(RegulatorData)
  Data_rownames=rownames(Data)
  
  # main loop
  # We want to end with the regulatory programs, hence we loop for Reassign, regulatory
  
  #STEP 1:  learning the regulatory program for each cluster
  ptm <- proc.time()
  regulatoryPrograms <- LINKER_LearnRegulatoryPrograms(Data,Clusters,RegulatorData,Lambda=Paramters$Lambda,alpha=Parameters$alpha,pmax=Parameters$pmax, mode=Parameters$mode, used_method=Parameters$used_method)
  ptm <- proc.time() - ptm
  #printf("Elapsed time is %f seconds\n",ptm[3])
  #print(rle(Clusters)$length)
  
  jj<-1
  while (jj < NrIterations){
    
    #STEP 2: reassigning genes based on closed match to new regulatory programs
    ptm <- proc.time()
    ReassignGenesToClusters <- LINKER_ReassignGenesToClusters(Data,RegulatorData,regulatoryPrograms$Beta,Clusters)
    ptm <- proc.time() - ptm
    #printf("Elapsed time at iteration %d is %f seconds\n",jj, ptm[3])
    jj<-jj+1
    
    NrReassignGenes = ReassignGenesToClusters$NrReassignGenes
    Clusters = ReassignGenesToClusters$Clusters
    #printf("Nr of reassignments is: %i \n\n",NrReassignGenes)
    
    #STEP 1:  learning the regulatory program for each cluster
    ptm <- proc.time()
    regulatoryPrograms <- LINKER_LearnRegulatoryPrograms(Data,Clusters,RegulatorData,Lambda=Paramters$Lambda,alpha=Parameters$alpha,pmax=Parameters$pmax, mode=Parameters$mode, used_method=Parameters$used_method)
    ptm <- proc.time() - ptm
    #printf("Elapsed time is %f seconds\n",ptm[3])
    #print(rle(Clusters)$length)
    
    # NrClusters = length(unique(Clusters))
    # sum = 0
    # for(i in 1:NrClusters){
    #   sum = sum + Matrix::nnzero(regulatoryPrograms$Beta[i,] )
    # }
    # avg = sum / NrClusters
    
    #printf("Average nr of regulators per module: %f \n",avg)
    
    #PreviousClusters = Clusters # using the clusters where the regulatory program was trained and not the last clusters
    
    
  }
  ptm1<- proc.time() - ptm1
  #printf("Elapsed time  is %f seconds\n\n",ptm1[3])
  
  # update results structure
  #ModuleMembership=as.matrix(PreviousClusters)
  ModuleMembership=as.matrix(Clusters)
  rownames(ModuleMembership)=rownames(Data)
  colnames(ModuleMembership)=c("ModuleNr")
  
  result <- list(NrModules = length(unique(Clusters)),RegulatoryPrograms = regulatoryPrograms$Beta,AllRegulators=rownames(RegulatorData),
                 AllGenes = rownames(Data),ModuleMembership = ModuleMembership)
  
  training_stats<-LINKER_EvaluateTestSet(result,Data,RegulatorData, used_method = LINKERinit$Parameters$used_method)
  
  result <- list(NrModules = length(unique(Clusters)),RegulatoryPrograms = regulatoryPrograms$Beta,AllRegulators=rownames(RegulatorData),
                 AllGenes = rownames(Data),ModuleMembership = ModuleMembership, trainingStats = training_stats)
  
  return(result)
}

LINKER_LearnRegulatoryPrograms<-function(Data,Clusters,RegulatorData,RegulatorSign,Lambda,alpha,pmax, mode, used_method="LINKER"){
  
  RegulatorData_rownames=rownames(RegulatorData)
  Data_rownames=rownames(Data)
  
  # stop has to be set because otherwise the algorithm continues until every
  # var is entered into the model
  #pmax = -10 # maximum nr of regulators that you want
  #trace = 0
  #NrFolds = 10
  NrClusters = length(unique(Clusters))
  NrGenes = nrow(Data)
  NrSamples = ncol(Data)
  #NrInterpolateSteps = 100
  
  y_all = mat.or.vec(NrClusters,NrSamples)
  ClusterIDs = unique(Clusters)
  ClusterIDs = sort(ClusterIDs, decreasing = FALSE)
  cnt <- 1:NrClusters
  
  # Compute the eigengen for each cluster
  ptm1 <- proc.time()
  #collapsed_rows<-collapseRows(Data,Clusters,rownames(Data),method='ME')
  #Data_collapsed<-collapsed_rows$datETcollapsed
  
  BetaY_all <- foreach(i=1:NrClusters,.combine=cbind,.init=list(list(),list(),list()),.packages = "glmnet") %dopar% {
    #for (i in 1:NrClusters){
    
    #y<-Data_collapsed[which(rownames(collapsed_rows$datETcollapsed)==as.character(ClusterIDs[i])),]
    
    CurrentClusterPositions = which(Clusters %in% ClusterIDs[i])
    nrGenesInClusters = length(CurrentClusterPositions)
    
    if(nrGenesInClusters>1)
    {
      if(used_method=="LINKER"){
        gene_svd<-svd(Data[CurrentClusterPositions,])
        y<-gene_svd$v[,1]
        if(cor(apply(Data[CurrentClusterPositions,],2,mean), y)<0){
          y<- -y
        }
        #y<-scale(y, center = FALSE)
      }
      else{
        y<-apply(Data[CurrentClusterPositions,],2,mean)
      }
    }
    else{
      y<-Data[CurrentClusterPositions,]
    }
    X = RegulatorData
    
    if(mode=="LASSO"){
      
      fit = cv.glmnet(t(X), y, alpha = alpha, pmax=10)
      #fit = cv.glmnet(t(X), y, alpha = alpha, pmax = pmax, lambda = grid)
      
      nonZeroLambdas <- fit$lambda[which(fit$nzero>0)]
      nonZeroCVMs <- fit$cvm[which(fit$nzero>0)]
      
      if(length(which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE)))==0){
        
        #for now: just print a warning, *although* this error WILL cause LINKER to crash in a few steps.
        warnMessage <- paste0("\nOn cluster ",i," there were no cv.glm results that gave non-zero coefficients.")
        warning(warnMessage);
        bestNonZeroLambda<-fit$lambda.min
      }
      else{
        bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE))]
      }
      #b_o = coef(fit,s = bestNonZeroLambda)
      b_o = coef(fit,s = fit$lambda.min)
      b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
    }
    else if(mode=="VBSR"){
      
      res<-vbsr(y,t(X),n_orderings = 15,family='normal')
      betas<-res$beta
      max_beta<-max(abs(betas))
      #betas[which(abs(betas)<0.1*max_beta)]<-0
      betas[res$pval > 0.05/(nrow(RegulatorData)*nrow(Data))]<-0
      b_opt<-betas
      b_o <- 0
      
    }
    
    else if(mode=="LM"){
      b_opt<-numeric(length = nrow(RegulatorData))
      for(i in 1:nrow(RegulatorData))
      {
        x<-X[i,]
        fit = lm(y~x)
        s<-summary(fit)
        if(s$coefficients[2,"Pr(>|t|)"]<0.05/(nrow(RegulatorData)*nrow(Data)))
        {
          b_opt[i]<-s$coefficients[2,1]
        }
        else{
          b_opt[i]<-0
        }
      }
      b_o <- 0
    }
    else{
      print("MODE NOT RECOGNIZED")
    }
    
    #y_all[i,] = y
    list(b_opt,y, b_o[1])
  }
  
  #ptm1<- proc.time() - ptm1
  #printf("Elapsed time is %f seconds\n",ptm1[3])
  
  tmpPos=NrClusters+1
  
  Beta <- do.call(cbind, BetaY_all[1,2:tmpPos])
  Beta = t(Beta);
  colnames(Beta)=RegulatorData_rownames
  rownames(Beta)=gsub('result.','Module_',rownames(Beta))
  
  y_all<-do.call(cbind, BetaY_all[2,2:tmpPos])
  y_all = t(y_all);
  rownames(y_all)=gsub('result.','Module_',rownames(y_all))
  
  intercept<-do.call(cbind, BetaY_all[3,2:tmpPos])
  intercept = t(intercept);
  intercept<-as.numeric(intercept)
  
  # calculating some statistics
  
  prediction<-(Beta %*% RegulatorData + intercept)
  error = y_all - prediction
  
  result <- list(Beta = Beta,error = error)
  return(result)
}

LINKER_ReassignGenesToClusters <- function(Data,RegulatorData,Beta,Clusters){
  
  MIN_NUM_GENES_PER_MODULE<-2
  
  RegulatorData_rownames=rownames(RegulatorData)
  Data_rownames=rownames(Data)
  
  NrGenes = nrow(Data)
  NrSamples = ncol(Data)
  NrReassignGenes = 0
  
  ##reassigning genes based on the Beta
  #getting the predictor data
  X = RegulatorData
  # creating the cluster "centroids"
  X1 = data.matrix(X)
  ModuleVectors = Beta %*% X1
  GeneNames = rownames(Data)
  
  #reassigning genes:
  #NewClusters = mat.or.vec(NrGenes,1)
  ptm1<- proc.time();
  nc <- foreach(i = 1:NrGenes, .combine = c) %dopar% {
    #for (i in 1:NrGenes){
    OldModule = Clusters[i]
    CurrentGeneVector = Data[i,,drop=FALSE]
    Correlations = abs(cor(t(CurrentGeneVector),t(ModuleVectors)))
    Correlations = (cor(t(CurrentGeneVector),t(ModuleVectors)))
    
    corr = data.matrix(Correlations,rownames.force = NA)
    #corr[is.na(corr)] <- -100000
    MaxCorrelation = max(corr,na.rm=TRUE)
    #print(MaxCorrelation)
    #MaxPosition = which(corr == max(corr,na.rm=TRUE))
    MaxPosition = which(signif(corr,digits=7) == signif(MaxCorrelation,digits=7))
    #print(MaxPosition)
    MaxPosition = MaxPosition[1] # this is new, to avoid two different reassignements
    
    if (MaxPosition != OldModule){
      NrReassignGenes = NrReassignGenes + 1
    }
    #NewClusters[i] = MaxPosition
    NewClusters = MaxPosition
    
  }
  
  # Remove cluster with too few genes. Avoids singularities. Could be solved imposing priors. future work
  for(i in unique(nc)){
    NrGenesInCluster<-sum(nc==i)
    if(NrGenesInCluster<MIN_NUM_GENES_PER_MODULE){
      # I need to reassign these genes
      genesInModule<-which(nc==i)
      print("One cluster removed!")
      #remove the cluster
      ModuleVectors[i,]<-0
      for(j in genesInModule){
        CurrentGeneVector = Data[j,,drop=FALSE]
        Correlations = abs(cor(t(CurrentGeneVector),t(ModuleVectors)))
        Correlations = (cor(t(CurrentGeneVector),t(ModuleVectors)))
        
        corr = data.matrix(Correlations,rownames.force = NA)
        MaxCorrelation = max(corr,na.rm=TRUE)
        MaxPosition = which(signif(corr,digits=7) == signif(MaxCorrelation,digits=7))
        MaxPosition = MaxPosition[1] # this is new, to avoid two different reassignements

        nc[j] = MaxPosition
        
      }
    }
  }
  
  ptm1<- proc.time() - ptm1;
  NrReassignGenes = length(which(nc!=Clusters));
  #print(NrReassignGenes)
  result <- list(NrReassignGenes = NrReassignGenes,Clusters = nc)
  return(result)
}

LINKER_EvaluateTestSet <- function(LINKERresults,MA_Data_TestSet,RegulatorData_TestSet, used_method="LINKER") {
  nrSamples = ncol(MA_Data_TestSet)
  RegulatorNames=rownames(RegulatorData_TestSet)
  
  #Iterating over the Modules
  stats = mat.or.vec(LINKERresults$NrModules,7)
  Rsquare = mat.or.vec(LINKERresults$NrModules,1)
  RsquareAjusted = mat.or.vec(LINKERresults$NrModules,1)
  modules <- list()
  
  for (i in 1:LINKERresults$NrModules){
    #check regulator presence
    currentRegulators = RegulatorNames[which(LINKERresults$RegulatoryPrograms[i,] != 0)]
    stats[i,1] = length(currentRegulators)
    
    #checking the presence of the clusters
    currentClusterGenes = LINKERresults$AllGenes[which(LINKERresults$ModuleMembership[,1] == i)]
    stats[i,2] = length(currentClusterGenes)
    
    #predict cluster expression in test set, always calculate but report
    #the totel percentage weight that is represented
    currentWeights = LINKERresults$RegulatoryPrograms[i,which(LINKERresults$RegulatoryPrograms[i,] != 0)]
    
    modules[[i]] = currentClusterGenes[currentClusterGenes %in% rownames(MA_Data_TestSet)]
    
    # drop=FALSE, this solves the problem when you have only one regulator, so the previous version is not needed.
    predictions = (t(RegulatorData_TestSet[currentRegulators,,drop=FALSE])) %*% (currentWeights) # need to make sure that the first argument remains a matrix.
    predictions = data.matrix(predictions)
      if (length(modules[[i]]) !=0) {
        if (length(currentClusterGenes)>1){
          cx <- MA_Data_TestSet[currentClusterGenes,]
          module_SVD = svd(cx)
          
          if(used_method=="LINKER"){
            outcome <- module_SVD$v[,1]
            if(cor(predictions,outcome)<0)
            {
              outcome <- -outcome
            }
            #outcome<-scale(outcome, center = FALSE)
          }
          else{
            outcome<-apply(cx,2,mean)
          }
          
          varEx<-module_SVD$d[1]^2/sum(module_SVD$d^2)
        } else {
          outcome = MA_Data_TestSet[currentClusterGenes,]
          varEx<-0
        }
        
        module_data<-MA_Data_TestSet[currentClusterGenes,]
        if(nrow(t(module_data))== 1){
          module_data<-t(module_data)
        }
        inmodule_corr<-abs(cor(t(module_data),outcome))
        meanInModuleCor<-mean(inmodule_corr)
        
        homogeneity<-abs(cor(t(module_data),t(module_data)))
        homogeneity<-(sum(homogeneity)-dim(homogeneity)[1])/((dim(homogeneity)[1]-1)*(dim(homogeneity)[1]-1))
        
        # using explained variance as metric, since mean square error is not
        # enough, no baseline interpretation possible
        
        SStot = sum((outcome-mean(outcome))^2)
        SSres = sum((predictions-outcome)^2)
        Rsquare = 1 - (SSres / SStot)
        RsquareAjusted = 1 - (1-Rsquare)*( (nrSamples-1)/(nrSamples - length(currentRegulators)))
        
        stats[i,3] = meanInModuleCor
        stats[i,4] = Rsquare
        stats[i,5] = RsquareAjusted
        stats[i,6]<-homogeneity
        stats[i,7]<- varEx
      } else {
        
      }
  }
  dimnames(stats) <- list(rownames(stats, do.NULL = FALSE, prefix = "Module_"), c("nrReg" , "nrGen", "MeanInModuleCorr", "Rsquare","RsquareAdjusted", "homogeneity", "condition"))

  return(stats)
}

LINKER_filterModuleBootstraps<-function(bootstrap_results,EvaluateTestSet)
{
  
  threashold<-0 # Let more clusters in
  field_to_filter<-3 #Use the R2 to filter rather than adj R2
  bootstrap_matrix<-sapply(bootstrap_results, function(x) x$ModuleMembership)
  rownames(bootstrap_matrix)<-bootstrap_results[[1]]$AllGenes
  
  filterField<-sapply(EvaluateTestSet, function(x) x[,field_to_filter])
  
  bad_modules<-sapply(filterField, function(x) which(x<threashold))
  NrGenes<-dim(bootstrap_matrix)[1]
  NrBootstraps<-dim(bootstrap_matrix)[2]
  for(i in 1:NrBootstraps)
  {
    idx<-which(bootstrap_matrix[, i] %in% bad_modules[[i]])
    bootstrap_matrix[idx, i]<-0
  }
  
  # idx<-1
  # genesToRemove<-numeric()
  # for(i in 1:NrGenes)
  # {
  #   if(sum(bootstrap_matrix[i,]==0)>0.3*NrBootstraps)
  #   {
  #     genesToRemove[idx]<-i
  #     idx<-idx+1
  #   }
  # }
  # if(idx>1)
  # {
  #   bootstrap_matrix<-bootstrap_matrix[-genesToRemove,]
  # }
  
  #print(idx)
  return(bootstrap_matrix)
}


LINKER_consensusClustering<-function(bootstrap_matrix)
{
  
  library(plyr)
  
  NrIter<-100
  
  consensus_clusters<-list()
  idx<-1
  
  NrGenes<-dim(bootstrap_matrix)[1]
  NrBootstraps<-dim(bootstrap_matrix)[2]
  consensusClust<-numeric(NrGenes)
  NrClusterVec<-apply(bootstrap_matrix,2,function(x)length(unique(x)))
  largest_bootstrap_idx<-which(NrClusterVec==max(NrClusterVec))[1]
  NrClusters<-max(NrClusterVec)
  K<-floor(NrGenes/NrClusters)
  
  gene_idxs<-1:NrGenes
  # for(i in 1:(NrClusters-1))
  # {
  #   gene_sampled<-sample(gene_idxs,K,replace = FALSE)
  #   consensusClust[gene_sampled]<-i
  #   gene_idxs<-gene_idxs[-which(gene_idxs %in% gene_sampled)]
  # }
  
  
  consensusClust<-bootstrap_matrix[,largest_bootstrap_idx]
  
  #consensusClust[gene_idxs]<-NrClusters
  cluster_rep<-matrix(0,nrow = NrClusters, ncol = NrBootstraps)
  
  for (iter in 1:NrIter)
  {
    
    # Compute Cluster Rep
    
    ClustersToRemove<-numeric(0)
    ctr<-1
    for (i in 1:NrClusters)
    {
      #get genes from cluster i
      cluster_genes<-which(consensusClust == i)
      if(length(cluster_genes) < 2)
      {
        ClustersToRemove[ctr]<-i
        ctr<-ctr+1
        next
      }
      cluster_matrix<-bootstrap_matrix[cluster_genes,]
      cluster_counts<-apply(cluster_matrix, 2, count)
      cluster_rep[i,]<-sapply(cluster_counts, function(x) x[which.max(x$freq),1])
    }
    if(length(ClustersToRemove)>1)
    {
      cluster_rep<-cluster_rep[-ClustersToRemove,]
      NrClusters<-nrow(cluster_rep)
    }
    
    #Assign genes to clusters
    for (gene_idx in 1:NrGenes)
    {
      min_dist<-NrBootstraps+1
      for (i in 1:NrClusters)
      {
        #compute distance with rep i
        tmp_dist<-length(which(cluster_rep[i,] != bootstrap_matrix[gene_idx,]))
        if(tmp_dist<min_dist){
          min_dist<-tmp_dist
          consensusClust[gene_idx]<-i
        }
      }
    }
    
  }
  
  # cluster_metric<-numeric(NrClusters)
  # cluster_sizes<-numeric(NrClusters)
  # for (i in 1:NrClusters)
  # {
  #   cluster_genes<-which(consensusClust == i)
  #   cluster_matrix<-bootstrap_matrix[cluster_genes,]
  #   hamm_dist<-apply(cluster_matrix, 1, function(x) length(which(cluster_rep[i,] != x)))
  #   cluster_metric[i]<-mean(hamm_dist)
  #   cluster_sizes[i]<-length(cluster_genes)
  #   
  #   if(cluster_metric[i] < 4)
  #   {
  #     bad_genes<-which(hamm_dist > 7)
  #     if(length(bad_genes)>0)
  #     {
  #       consensus_clusters[[idx]]<-rownames(cluster_matrix[-bad_genes,])
  #     }
  #     else
  #     {
  #       consensus_clusters[[idx]]<-rownames(cluster_matrix)
  #     }
  #     
  #     idx<-idx+1
  #   }
  # }
  
  return(consensusClust)
}

LINKER_consensusClustering_IPC<-function(bootstrap_matrix)
{
  
  library(plyr)
  
  NrIter<-100
  
  consensus_clusters<-list()
  idx<-1
  
  NrGenes<-dim(bootstrap_matrix)[1]
  NrBootstraps<-dim(bootstrap_matrix)[2]
  consensusClust<-numeric(NrGenes)
  newConsensusClust<-numeric(NrGenes)
  NrClusterVec<-apply(bootstrap_matrix,2,function(x)length(unique(x)))
  largest_bootstrap_idx<-which(NrClusterVec==max(NrClusterVec))[1]

  rnd_bootstrap_idx<-sample(NrBootstraps, 1)
  consensusClust<-bootstrap_matrix[,rnd_bootstrap_idx]
  
  clustID<-unique(consensusClust)
  NrClusters<-length(clustID)
  
  gene_idxs<-1:NrGenes
  
  similarity_matrix<-matrix(nrow = NrGenes, ncol = NrGenes)
  for(i in 1:NrGenes){
    for(j in 1:NrGenes){
      similarity_matrix[i,j]<-sum(bootstrap_matrix[i,]==bootstrap_matrix[j,])
    }
  }
  
  for (iter in 1:NrIter)
  {
    #Assign genes to clusters
    for (gene_idx in 1:NrGenes)
    {
      # assign gene to most similar cluster
      max_dist<-0
      temp_dist<-0
      for (i in 1:NrClusters)
      {
        #compute distance with cluster i
        cluster_genes<-which(consensusClust == clustID[i])
        if(length(cluster_genes)<2){
          #print("foooooo")
        }
        tmp_dist<-sum(apply(as.matrix(cluster_genes), 1, function(x)similarity_matrix[x,gene_idx]))/length(cluster_genes)
        if(tmp_dist>max_dist){
          max_dist<-tmp_dist
          newConsensusClust[gene_idx]<-i
        }
      }
    }
    if(sum(consensusClust==newConsensusClust)==length(consensusClust)){
      break
    }
    consensusClust<-newConsensusClust
    clustID<-unique(consensusClust)
    NrClusters<-length(clustID)
    
  }
  
 
  
  return(consensusClust)
}

LINKER_compute_lincRegulators<-function(Data, RegulatorData, Clusters, Parameters)
{
  
  # valid_genes<-numeric()
  # NrModules<-length(consensus_modules)
  # for(i in 1:length(consensus_modules))
  # {
  #   valid_genes<-c(valid_genes,which(rownames(Data) %in% consensus_modules[[i]]))
  # }
  # Data<-Data[valid_genes,]
  # 
  # Clusters<-numeric(length=nrow(Data))
  # names(Clusters) <- rownames(Data)
  # for(i in 1:length(consensus_modules))
  # {
  #   module_genes<-which(rownames(Data) %in% consensus_modules[[i]])
  #   Clusters[module_genes]<-i
  # }
  
  regulatoryPrograms_results <- LINKER_LearnRegulatoryPrograms(Data,Clusters,RegulatorData,Lambda=Paramters$Lambda,alpha=Parameters$alpha,pmax=Parameters$pmax, mode=Parameters$mode, used_method=Parameters$used_method)
  
  
  Results<-list(NrModules=length(unique(Clusters)), RegulatoryPrograms=regulatoryPrograms_results$Beta, ModuleMembership=as.matrix(Clusters), AllGenes=rownames(Data))
  
  return(Results)
}


LINKER_createBootstrapGraph<-function(bootstrap_results, bootstrap_modules)
{
  
  over_enrichment_pvalues<-numeric()
  p_idx<-1
  total_genes<-10000
  
  Nr_bootstraps<-length(bootstrap_results)
  
  for(i in 1: Nr_bootstraps){
    
    Nr_modules_source<-bootstrap_results[[i]]$NrModules
    
    ## compare the old bootsrap modules with the modules on the rest of bootstraps
    for(j in 1:Nr_modules_source){
      
      for(ii in 1: Nr_bootstraps){
        
        Nr_modules_sink<-bootstrap_results[[ii]]$NrModules
        
        for(jj in 1:Nr_modules_sink){
          
          white_balls<-length(bootstrap_modules[[i]][[j]])
          black_balls<-total_genes - white_balls
          
          drawn<-length(bootstrap_modules[[ii]][[jj]])
          drawn_whites<-length(which(bootstrap_modules[[i]][[j]] %in% bootstrap_modules[[ii]][[jj]]))
          
          over_enrichment_pvalues[p_idx]<- phyper(drawn_whites-1, white_balls, black_balls, drawn, lower.tail = FALSE, log.p = FALSE)
          p_idx<-p_idx+1
        }
      }
    }
  }
  
  over_adjp<-p.adjust(over_enrichment_pvalues,'holm')
  
  over_adjp_Matrix<- matrix(over_adjp, nrow = sqrt(length(over_adjp)), byrow = TRUE)
  
  enriched_idx<-which(over_adjp_Matrix<1e-10,arr.ind=TRUE)
  
  adj_sparseMatrix_graph<-sparseMatrix(enriched_idx[,1], enriched_idx[,2])
  
  bootstrap_graph<-graph_from_adjacency_matrix(as.matrix(adj_sparseMatrix_graph), mode = "directed",  diag = FALSE)
  
  bootstrap_comm<-cluster_walktrap(bootstrap_graph, weights = E(bootstrap_graph)$weight, steps = 4,
                                   merges = TRUE, modularity = TRUE, membership = TRUE)
  
  return(list(bootstrap_graph, bootstrap_comm))
}

LINKER_createEnrichedGraph<-function(enriched_modules)
{
  
  over_enrichment_pvalues<-numeric()
  p_idx<-1
  total_genes<-10000
  

    
    Nr_modules<-length(enriched_modules)
    
    ## compare the old bootsrap modules with the modules on the rest of bootstraps
  for(j in 1:Nr_modules){
      
      for(jj in 1:Nr_modules){
        
        white_balls<-length(enriched_modules[[j]]$protein_coding_genes)
        black_balls<-total_genes - white_balls
        
        drawn<-length(enriched_modules[[jj]]$protein_coding_genes)
        drawn_whites<-length(which(enriched_modules[[jj]]$protein_coding_genes %in% enriched_modules[[j]]$protein_coding_genes))
        
        over_enrichment_pvalues[p_idx]<- phyper(drawn_whites-1, white_balls, black_balls, drawn, lower.tail = FALSE, log.p = FALSE)
        p_idx<-p_idx+1
      }
  }
  
  over_adjp<-p.adjust(over_enrichment_pvalues,'holm')
  
  over_adjp_Matrix<- matrix(over_adjp, nrow = sqrt(length(over_adjp)), byrow = TRUE)
  
  enriched_idx<-which(over_adjp_Matrix<0.05,arr.ind=TRUE)
  
  adj_sparseMatrix_graph<-sparseMatrix(enriched_idx[,1], enriched_idx[,2])
  
  bootstrap_graph<-graph_from_adjacency_matrix(as.matrix(adj_sparseMatrix_graph), mode = "undirected",  diag = FALSE)
  
  bootstrap_comm<-cluster_walktrap(bootstrap_graph, weights = E(bootstrap_graph)$weight, steps = 4,
                                   merges = TRUE, modularity = TRUE, membership = TRUE)
  
  return(list(bootstrap_graph, bootstrap_comm))
}

LINKER_computeBootstrapModuleEnrichment<-function(old_modules, bootstrap_modules)
{
  
  
  
  #### Compare the new results to the old modules
  
  
  return(enriched_modules)
  
}

LINKER_run<-function(lognorm_est_counts, protein_filtered_idx, lincs_filtered_idx, Gene_set_Collections,
                     link_mode=c("LASSO", "VBSR"),
                     module_rep="MEAN",
                     NrModules=100, 
                     corrClustNrIter=100,
                     Nr_bootstraps=10,
                     FDR=0.05,
                     module_summary=0)
{
  res<-list()
  modules<-list()

  for(i in 1:length(link_mode)){
    res[[ link_mode[i] ]]<-run_linker(lognorm_est_counts, protein_filtered_idx,  lincs_filtered_idx, NrModules, module_summary,
                                      mode=link_mode[i], used_method=module_rep, corrClustNrIter,Nr_bootstraps)
    modules[[ link_mode[i] ]]<-filter_enriched_modules(Gene_set_Collections,res[[ link_mode[i] ]],FDR)
  }
  
  graphs<-LINKER_compute_graphs_from_modules(modules, lognorm_est_counts)
  

}
