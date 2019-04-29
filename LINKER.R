LINKER_run<-function(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx, Gene_set_Collections,
                     link_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                     graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                     module_rep="MEAN",
                     NrModules=100, 
                     corrClustNrIter=100,
                     Nr_bootstraps=10,
                     FDR=0.05,
                     NrCores=30)
  
{
  res<-list()
  modules<-list()
  graphs<-list()
  for(i in 1:length(link_mode)){
    res[[ link_mode[i] ]]<-LINKER_runPhase1(lognorm_est_counts, target_filtered_idx,  regulator_filtered_idx, 
                                      NrModules,NrCores=NrCores,
                                      mode=link_mode[i], used_method=module_rep, 
                                      corrClustNrIter=corrClustNrIter,Nr_bootstraps=Nr_bootstraps)
    
    modules[[ link_mode[i] ]]<-LINKER_extract_modules(res[[ link_mode[i] ]])
    print(paste0("Link mode ",link_mode[i]," completed!"))  
    
    graphs[[ link_mode[i] ]]<-list()
    for(j in 1:length(graph_mode)){
      graphs[[ link_mode[i] ]][[ graph_mode[j] ]] <- LINKER_compute_modules_graph(modules[[ link_mode[i] ]], lognorm_est_counts, mode=graph_mode[j])
      print(paste0("Graphs for (",link_mode[i],",",graph_mode[j], ") computed!"))    
    }
    
  }
  
  GEAs<-LINKER_compute_graph_enrichment_geneSets_graph_list(Gene_set_Collections,graphs,FDR=FDR,BC=1, NrCores = NrCores)
  
  return(list(raw_results=res,modules=modules,graphs=graphs, GEAs=GEAs))
  
  
}

############# PHASE 1 functions ###################

LINKER_runPhase1<-function(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx, NrModules, 
                     Lambda=0.0001, alpha=1-1e-06,
                     pmax=10, mode="LASSO", used_method="LINKER",
                     NrCores=30, corrClustNrIter=21,
                     Nr_bootstraps=1)
{
  
  # Creating the parameters structure
  Parameters <- list(Lambda=Lambda,pmax=pmax,alpha=alpha, mode=mode, used_method=used_method)
  sample_size<-dim(lognorm_est_counts)[2]
  #test_size<-round(0.2*sample_size)
  train_size<-round(0.8*sample_size)
  EvaluateTestSet<-list()
  bootstrap_modules<-list()
  bootstrap_results<-list()
  #test_samples<-sample(1:sample_size, test_size, replace=F)
  #train_val_samples<-setdiff(1:sample_size, test_samples)
  
  #Regulator_data_test = t(scale(t(lognorm_est_counts[regulator_filtered_idx,test_samples])))
  #MA_matrix_Var_test = t(scale(t(lognorm_est_counts[target_filtered_idx,test_samples])))
  
  #Regulator_data_train_val = t(scale(t(lognorm_est_counts[regulator_filtered_idx,train_val_samples])))
  #MA_matrix_Var_train_val = t(scale(t(lognorm_est_counts[target_filtered_idx,train_val_samples])))
  
  for(boost_idx in 1:Nr_bootstraps)
  {
    
    train_samples<-sample(1:sample_size, train_size, replace=F)
    validation_samples<-setdiff(1:sample_size, train_samples)
    
    Regulator_data_train = t(scale(t(lognorm_est_counts[regulator_filtered_idx,train_samples])))    
    Regulator_data_validation = t(scale(t(lognorm_est_counts[regulator_filtered_idx,validation_samples])))
    
    MA_matrix_Var_train = t(scale(t(lognorm_est_counts[target_filtered_idx,train_samples])))
    MA_matrix_Var_validation = t(scale(t(lognorm_est_counts[target_filtered_idx,validation_samples])))
    
    LINKERinit<-LINKER_init(MA_matrix_Var = MA_matrix_Var_train, RegulatorData = Regulator_data_train, NrModules = NrModules, NrCores=NrCores, corrClustNrIter=corrClustNrIter, Parameters = Parameters )
    
    
    tmp<-LINKER_corrClust(LINKERinit)
    bootstrap_results[[boost_idx]]<-tmp
    
    EvaluateTestSet[[boost_idx]]<-LINKER_EvaluateTestSet(bootstrap_results[[boost_idx]],MA_matrix_Var_validation,Regulator_data_validation, used_method=Parameters$used_method)
    
    printf("Bootstrap %d, NrModules %d:\n", boost_idx, bootstrap_results[[boost_idx]]$NrModules)
    
    #print(EvaluateTestSet[[boost_idx]])
    
    print(apply(EvaluateTestSet[[boost_idx]], 2, mean))
  }
  
  
  return(list(bootstrapResults=bootstrap_results,bootstrapTestStats=EvaluateTestSet))
  
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
  
  BetaY_all <- foreach(i=1:NrClusters,.combine=cbind,.init=list(list(),list(),list()),.packages = c("vbsr","glmnet","igraph")) %dopar% {
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
    
    if(mode=="LASSOmin"){
      
      fit = cv.glmnet(t(X), y, alpha = alpha)
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
    else if(mode=="LASSO1se"){
      
      fit = cv.glmnet(t(X), y, alpha = alpha)
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
      b_o = coef(fit,s = fit$lambda.1se)
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
      #print("One cluster removed!")
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

LINKER_extract_modules<-function(results){
  
  modules<-list()
  
  enriched_idx<-1
  NrBootstraps<-length(results$bootstrapResult)
  for(idx_bootstrap in 1:NrBootstraps){
    
    NrModules<-results$bootstrapResult[[idx_bootstrap]]$NrModules
    boot_results<-results$bootstrapResult[[idx_bootstrap]]
    boot_idx<-sort(unique(boot_results$ModuleMembership[,]))
    
    for(Module_number in 1:NrModules){

      
      Module_target_genes_full_name<-boot_results$AllGenes[which(boot_results$ModuleMembership[,]==boot_idx[Module_number])]
      Module_target_gene_list<-sapply(Module_target_genes_full_name, function(x) strsplit(x, "\\|"))
      if(length(Module_target_gene_list[[1]])==1)
      {
        #NOT FULL GENECODE NAME, ONLY ONE NAME PER GENE!
        Module_target_genes<-Module_target_genes_full_name
      }
      else{
        #FULL GENECODE ANNOTATION!
        Module_target_genes<-sapply(Module_target_gene_list, function(x) x[[6]])
        Module_target_genes<-unname(Module_target_genes)
      }
      
      
      Modules_regulators_full_name<-names(which(boot_results$RegulatoryPrograms[Module_number,]!=0))
      if(length(Modules_regulators_full_name)==0){
        next;
      }
      Modules_regulators_list<-sapply(Modules_regulators_full_name, function(x) strsplit(x, "\\|"))
      if(length(Modules_regulators_list[[1]])==1)
      {
        #NOT FULL GENECODE NAME, ONLY ONE NAME PER GENE!
        Modules_regulators<-Modules_regulators_full_name
      }
      else{
        #FULL GENECODE ANNOTATION!
        Modules_regulators<-sapply(Modules_regulators_list, function(x) x[[6]])
        Modules_regulators<-unname(Modules_regulators)
      }
      
      
      
      modules[[enriched_idx]]<-list(
        target_genes=Module_target_genes_full_name, 
        regulators=Modules_regulators_full_name, 
        regulatory_program=boot_results$RegulatoryPrograms[Module_number,],
        training_stats=boot_results$trainingStats[Module_number,],
        test_stats=results$bootstrapTestStats[[idx_bootstrap]][Module_number],
        assigned_genes=which(boot_results$ModuleMembership[,]==Module_number),
        bootstrap_idx=idx_bootstrap
      )
      enriched_idx<-enriched_idx+1
    }
  }
  return(modules)
  
}

######### Phase2 Functions ##########

LINKER_compute_modules_graph<-function(modules, Data, mode="VBSR",alpha=1-1e-06)
{
  
  bp_g<-list()
  i<-1
  
  bp_g<-foreach(mod_idx=1:length(modules), .packages = c("vbsr","glmnet","igraph"))%dopar%
    #for(mod_idx in 1:length(modules))
  {
    targetgenes<-unlist(modules[[mod_idx]]$target_genes)
    regulators<-unlist(modules[[mod_idx]]$regulators)
    X<-Data[regulators,]
    
    #We need to handle the special case where only one regulator regulates a module/community
    if(length(regulators)<2){
      
      non_zero_beta<-modules[[mod_idx]]$regulatory_program[which(modules[[mod_idx]]$regulatory_program != 0)]
      if(length(non_zero_beta) != 1){
        warning("NON_ZERO_BETA != 1")
      }
      
      driverMat<-matrix(data = non_zero_beta, nrow = length(targetgenes), ncol = length(regulators))
      
      rownames(driverMat)<-targetgenes
      colnames(driverMat)<-regulators
    }
    else{
      
      driverMat<-matrix(data = NA, nrow = length(targetgenes), ncol = length(regulators))
      
      for(idx_gene in 1:length(targetgenes))
      {
        y<-Data[targetgenes[idx_gene],]
        
        if(mode=="VBSR")
        {
          res<-vbsr(y,t(X),n_orderings = 15,family='normal')
          betas<-res$beta
          betas[res$pval > 0.05/(length(regulators)*length(targetgenes))]<-0
          driverMat[idx_gene,]<-betas
        }
        else if(mode=="LASSOmin")
        {
          fit = cv.glmnet(t(X), y, alpha = alpha)
          
          b_o = coef(fit,s = fit$lambda.min)
          b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
          driverMat[idx_gene,]<-b_opt 
        }
        else if(mode=="LASSO1se")
        {
          fit = cv.glmnet(t(X), y, alpha = alpha)
          
          b_o = coef(fit,s = fit$lambda.1se)
          b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
          driverMat[idx_gene,]<-b_opt 
        }
        else if(mode=="LM")
        {
          for(idx_regs in 1:length(regulators))
          {
            x<-t(X)[,idx_regs]
            fit = lm(y~x)
            s<-summary(fit)
            driverMat[idx_gene,idx_regs]<-s$coefficients[2,"Pr(>|t|)"]<0.05/(length(targetgenes)*length(regulators))
          }
        }
        else
        {
          warning("MODE NOT RECOGNIZED")
        }
        
      }
      
      rownames(driverMat)<-targetgenes
      colnames(driverMat)<-regulators
      
      regulated_genes<-which(rowSums(abs(driverMat))!=0)
      regulatory_genes<-which(colSums(abs(driverMat))!=0)
      
      # We need to treat the special cases independently
      if(length(regulated_genes)<2){
        driverMat<-driverMat[regulated_genes,]
        driverMat<-driverMat[regulatory_genes]
      }
      else if(length(regulatory_genes)<2){
        driverMat<-driverMat[,regulatory_genes]
        driverMat<-driverMat[regulated_genes]
      }
      else{
        driverMat<-driverMat[,regulatory_genes]
        driverMat<-driverMat[regulated_genes,]
      }
      
    }
    
    #bp_g[[i]]<-
    graph_from_incidence_matrix(driverMat)
    #i<-i+1
  }
  
  return(bp_g)
}

############## Gene enrichment functions #####################

#WARNING: The functions expect the gene annotation either using the full Genecode annotation or just the Gene name

LINKER_compute_graph_enrichment_geneSets_graph_list<-function(pathway_genes,g,FDR=0.05,BC=1, NrCores=1)
{
  
  GEA<-list()
  
  for(i in 1:length(g))
  {
    link_mode<-names(g)[i]
    GEA[[ link_mode ]]<-list()
    
    for(j in 1:length(g[[link_mode]]))
    {
      graph_mode<-names(g[[link_mode]])[j]
      GEA[[ link_mode ]][[graph_mode]]<-LINKER_compute_graph_list_enrichment_geneSets(pathway_genes,g[[link_mode]][[graph_mode]], BC=BC, NrCores=NrCores)
      print(paste0("GEA for (",link_mode,",",graph_mode, ") computed!"))      
    }
  }
  
  return(GEA)
}

LINKER_compute_graph_list_enrichment_geneSets<-function(pathway_genes,g_list,FDR=0.05,BC=1, NrCores=1)
{
  GEA<-list()
  
  Num_regs<-sapply(g_list, function(x) sum(V(x)$type==1))
  #BC<-mean(Num_regs)*BC
  # BIOCARTA
  Gene_set_Collections<-pathway_genes[1]
  GEA$BIOCARTA<-LINKER_compute_reg_enrichment_from_graph_list(g_list, Gene_set_Collections,FDR=FDR, BC=BC,NrCores=NrCores)
  # KEGG
  Gene_set_Collections<-pathway_genes[2]
  GEA$KEGG<-LINKER_compute_reg_enrichment_from_graph_list(g_list, Gene_set_Collections,FDR=FDR, BC=BC,NrCores=NrCores)
  # REACTOME
  Gene_set_Collections<-pathway_genes[3]
  GEA$REACTOME<-LINKER_compute_reg_enrichment_from_graph_list(g_list, Gene_set_Collections,FDR=FDR, BC=BC,NrCores=NrCores)
  # GENESIGDB
  Gene_set_Collections<-pathway_genes[4]
  GEA$GENESIGDB<-LINKER_compute_reg_enrichment_from_graph_list(g_list, Gene_set_Collections,FDR=FDR, BC=BC,NrCores=NrCores)
  # ALL
  #Gene_set_Collections<-pathway_genes[c(3,4,5,12)]
  #GEA$ALL<-LINKER_compute_reg_enrichment_from_graph_list(g_list, Gene_set_Collections,FDR=FDR, BC=BC)
  
  return(GEA)
}

LINKER_compute_reg_enrichment_from_graph_list<-function(g_list, Gene_set_Collections,FDR=0.05, BC=1,NrCores=1)
{
  
  registerDoParallel(NrCores)
  path_regs<-list()
  for(i in 1:length(g_list))
  {
    g<-g_list[[i]]
    reg_neighbors<-sapply(V(g)[V(g)$type==1], function(x) neighbors(g,x))
    GEA<-mclapply(reg_neighbors,function(x) LINKER_module_gene_enrichment(names(x), Gene_set_Collections, FDR, BC))
    path_regs[[i]]<-unlist(lapply(GEA,function(x) unlist(x)))
  }
  
  #path_regs<-unlist(path_regs)
  #unique_paths_reg<-unique(names(path_regs))
  #GEA_per_reg<-path_regs[unique_paths_reg]
  
  return(path_regs)
  
}

## This function computes the enrichment analysis from a set of genes
LINKER_module_gene_enrichment <-function(Module, pathway_genes, FDR=0.05, BC=1, total_genes=10000)
{
  
  i<-1
  #under_enrichment_pvalues<-numeric()
  over_enrichment_pvalues<-numeric()
  for(collection_idx in 1:length(pathway_genes)){
    
    for(pathway_idx in 1:length(pathway_genes[[collection_idx]])){
      
      white_balls<-length(pathway_genes[[collection_idx]][[pathway_idx]])
      black_balls<-total_genes - white_balls
      
      drawn<-length(Module)
      drawn_whites<-length(which(pathway_genes[[collection_idx]][[pathway_idx]] %in% Module))
      
      #under_enrichment_pvalues[i]<- phyper(drawn_whites, white_balls, black_balls, drawn, lower.tail = TRUE, log.p = FALSE)
      over_enrichment_pvalues[i]<- phyper(drawn_whites-1, white_balls, black_balls, drawn, lower.tail = FALSE, log.p = FALSE)
      i<-i+1
    }
  }
  
  if(BC<0){
    #under_adjp<-under_adjp*(-BC)
    over_adjp<-over_adjp*(-BC)
  }
  
  else{
    #under_adjp<-p.adjust(under_enrichment_pvalues,'fdr')
    over_adjp<-p.adjust(over_enrichment_pvalues,'fdr')
    #under_adjp<-under_adjp*BC
    over_adjp<-over_adjp*BC
  }
  
  
  
  i<-1
  GEA_over<-list()
  #GEA_under<-list()
  for(collection_idx in 1:length(pathway_genes))  
  {
    j<-1
    jj<-1
    #GEA_under[[collection_idx]]<-numeric()
    GEA_over[[collection_idx]]<-numeric()
    gensetNames_over<-character()
    #gensetNames_under<-character()
    for(pathway_idx in 1:length(pathway_genes[[collection_idx]]))
    {
      if(over_adjp[i]<FDR){
        GEA_over[[collection_idx]][j]<-over_adjp[i]
        gensetNames_over[j]<-names(pathway_genes[[collection_idx]])[pathway_idx]
        j<-j+1
      }
      
      # if(under_adjp[i]<FDR){
      #   GEA_under[[collection_idx]][jj]<-over_adjp[i]
      #   gensetNames_under[jj]<-names(pathway_genes[[collection_idx]])[pathway_idx]
      #   jj<-jj+1
      # }
      i<-i+1
    }
    names(GEA_over[[collection_idx]])<-gensetNames_over
    #names(GEA_under[[collection_idx]])<-gensetNames_under
  }
  #return(list(GEA_under, GEA_over))
  return(GEA_over)
}
##