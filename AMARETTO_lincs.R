AMARETTO_linc_init <- function(MA_matrix_Var, RegulatorData, NrModules,
                         PvalueThreshold=0.001,RsquareThreshold=0.1,pmax=10,NrCores=1,OneRunStop=0) {
  
  # Fixing default paramters
  AutoRegulation=2
  Lambda2=0.0001
  Lambda2=0.100
  alpha=1-1e-06
  
  # Creating the parameters structure
  Parameters <- list(AutoRegulation=AutoRegulation,OneRunStop=OneRunStop,Lambda2=Lambda2,Mode='larsen',pmax=pmax,alpha=alpha)
  
  if (nrow(MA_matrix_Var)>NrModules){
    #KmeansResults=kmeans(MA_matrix_Var,NrModules,iter.max=100)
    
    #Random initailization of the centroids
    rnd_cent<-runif(NrModules, min = 1, max = nrow(MA_matrix_Var))
    #rnd_cent<-1:NrModules
    ModuleVectors<-MA_matrix_Var[rnd_cent,]

    Data<-MA_matrix_Var
    Clusters<-numeric()
    for (i in 1:nrow(MA_matrix_Var)){
      CurrentGeneVector = Data[i,,drop=FALSE]
      Correlations = cor(t(CurrentGeneVector),t(ModuleVectors))
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
  } else {
    stop("The number of modules is too large compared to the total number of genes.")
  }
  #ModuleMembership=as.numeric(KmeansResults$cluster)
  ModuleMembership<-as.numeric(Clusters)
  names(ModuleMembership) <- rownames(MA_matrix_Var)
  
  return(list(MA_matrix_Var=MA_matrix_Var,RegulatorData=RegulatorData,ModuleMembership=ModuleMembership,Parameters=Parameters,NrCores=NrCores))

}
AMARETTO_linc_run <- function(AMARETTOinit)
{                              
  if (nrow(AMARETTOinit$RegulatorData)==1){
    stop('For cancer ',CancerSite,' only one driver is detected. AMARETTO cannot be run with less than two drivers.\n')
  }
  cat('Running AMARETTO on',length(rownames(AMARETTOinit$MA_matrix_Var)),'genes and',length(colnames(AMARETTOinit$MA_matrix_Var)),'samples.\n')
  cat('\tStopping if less then',0.01*length(rownames(AMARETTOinit$MA_matrix_Var)),'genes reassigned.\n')
  result=AMARETTO_LarsenBased_linc(AMARETTOinit$MA_matrix_Var,AMARETTOinit$ModuleMembership,AMARETTOinit$RegulatorData,AMARETTOinit$Parameters,AMARETTOinit$NrCores)
  
  #result$ModuleData=AMARETTO_CreateModuleData(AMARETTOinit,result)
  #result$RegulatoryProgramData=AMARETTO_CreateRegulatorPrograms(AMARETTOinit,result)
  
  return(result)
}

AMARETTO_LarsenBased_linc <- function(Data,Clusters,RegulatorData,Parameters,NrCores){
  # this will register nr of cores/threads, keep this here so the user can decide how many cores based on their hardware.
  registerDoParallel(cores=NrCores)
  ptm1 <- proc.time()
  
  NrIterations<-5

  
  RegulatorData_rownames=rownames(RegulatorData)
  Data_rownames=rownames(Data)
  
  AutoRegulation = Parameters$AutoRegulation
  RegulatorSign=array(0,length(RegulatorData_rownames))
  Lambda = Parameters$Lambda2
  OneRunStop = Parameters$OneRunStop
  if (AutoRegulation == 1){
    cat('\tAutoregulation is turned ON.\n')
  } else if (AutoRegulation == 2){
    cat('\tAutoregulation is turned ON.\n')
  } else {
    cat('\tAutoregulation is turned OFF.\n')
  }
  
  # main loop
  NrReassignGenes = length(Data_rownames)
  jj<-1
  #while (NrReassignGenes > 0.01*length(Data_rownames)){
  while (jj < NrIterations){
    
    #STEP 1:  learning the regulatory program for each cluster
    ptm <- proc.time()
    print(length(unique(Clusters)))
    switch(Parameters$Mode,
           larsen={
             regulatoryPrograms <- AMARETTO_LearnRegulatoryProgramsLarsen(Data,Clusters,RegulatorData,RegulatorSign,Lambda,AutoRegulation,alpha=Parameters$alpha,pmax=Parameters$pmax)
           }
    )
    ptm <- proc.time() - ptm
    printf("Elapsed time is %f seconds\n",ptm[3])
    
    NrClusters = length(unique(Clusters))
    sum = 0
    for(i in 1:NrClusters){
      sum = sum + Matrix::nnzero(regulatoryPrograms$Beta[i,] )
    }
    avg = sum / NrClusters
    
    printf("Average nr of regulators per module: %f \n",avg)
    
    PreviousClusters = Clusters # using the clusters where the regulatory program was trained and not the last clusters
    if (OneRunStop == 1){ break } 		# running only one iteration of optimization, useful in large comparisons
    
    #STEP 2: reassigning genes based on closed match to new regulatory programs
    ptm <- proc.time()
    ReassignGenesToClusters <- AMARETTO_ReassignGenesToClusters(Data,RegulatorData,regulatoryPrograms$Beta,Clusters,AutoRegulation)
    ptm <- proc.time() - ptm
    printf("Elapsed time at iteration %d is %f seconds\n",jj, ptm[3])
    jj<-jj+1
    
    NrReassignGenes = ReassignGenesToClusters$NrReassignGenes
    Clusters = ReassignGenesToClusters$Clusters
    printf("Nr of reassignments is: %i \n\n",NrReassignGenes)
  }
  ptm1<- proc.time() - ptm1
  printf("Elapsed time  is %f seconds\n\n",ptm1[3])
  
  # update results structure
  ModuleMembership=as.matrix(PreviousClusters)
  rownames(ModuleMembership)=rownames(Data)
  colnames(ModuleMembership)=c("ModuleNr")
  
  result <- list(NrModules = length(unique(Clusters)),RegulatoryPrograms = regulatoryPrograms$Beta,AllRegulators=rownames(RegulatorData),
                 AllGenes = rownames(Data),ModuleMembership = ModuleMembership,AutoRegulationReport=regulatoryPrograms$AutoRegulationReport)
  
  return(result)
}

AMARETTO_LearnRegulatoryProgramsLarsen_foo<-function(Data,Clusters,RegulatorData,RegulatorSign,Lambda,AutoRegulation,alpha,pmax){
  
  RegulatorData_rownames=rownames(RegulatorData)
  Data_rownames=rownames(Data)
  
  # stop has to be set because otherwise the algorithm continues until every
  # var is entered into the model
  #pmax = -10 # maximum nr of regulators that you want
  trace = 0
  NrFolds = 10
  NrClusters = length(unique(Clusters))
  NrGenes = nrow(Data)
  NrSamples = ncol(Data)
  NrInterpolateSteps = 100
  
  # autoregulation yes or no?
  if (AutoRegulation >= 1){
    #Beta = mat.or.vec((NrClusters),length(RegulatorData_rownames))
  } else if (AutoRegulation == 0) {
    BetaSpecial = list(NrClusters,1)
    RegulatorPositions = list(NrClusters,1)
  }
  #AutoRegulationReport = mat.or.vec(NrClusters,3)
  
  y_all = mat.or.vec(NrClusters,NrSamples)
  y_std = numeric(NrClusters)
  
  ClusterIDs = unique(Clusters)
  ClusterIDs = sort(ClusterIDs, decreasing = FALSE)
  cnt <- 1:NrClusters
  
  ptm1 <- proc.time()
  BetaY_all <- foreach(i=1:NrClusters,.combine=cbind,.init=list(list(),list(),list()),.packages = "glmnet") %dopar% {
    #for (i in 1:NrClusters){
    
    if (length(which(Clusters == ClusterIDs[i]))>1) {
      y = apply((Data[which(Clusters == ClusterIDs[i]),]),2,mean)
    } else {
      y = Data[which(Clusters == ClusterIDs[i]),]
    }
                 
    CurrentClusterPositions = which(Clusters %in% ClusterIDs[i])
    nrGenesInClusters = length(CurrentClusterPositions)
    
    if (AutoRegulation >= 1){
      X = RegulatorData
    } else if (AutoRegulation == 0){
      X = RegulatorData[setdiff(RegulatorData_rownames,Data_rownames[CurrentClusterPositions]),]
    }
    
    fit = cv.glmnet(t(X), y,alpha = alpha, pmax = pmax)
    
    nonZeroLambdas <- fit$lambda[which(fit$nzero>0)]
    nonZeroCVMs <- fit$cvm[which(fit$nzero>0)]
    
    if(length(which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE)))==0){
      
      #for now: just print a warning, *although* this error WILL cause Amaretto to crash in a few steps.
      warnMessage <- paste0("\nOn cluster ",i," there were no cv.glm results that gave non-zero coefficients.")
      warning(warnMessage);
      
    }
    
    bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE))]
    b_o = coef(fit,s = bestNonZeroLambda)
    b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
    
    if (AutoRegulation == 2){ # autoregulation is allowed, if the regulator is stable also after removing it from the its regulated cluster.
      
      CurrentUsedRegulators = RegulatorData_rownames[which(b_opt!=0, arr.ind = T)]
      CurrentClusterMembers = Data_rownames[CurrentClusterPositions]
      nrIterations = 0 # 0 means no overlap initially
      while (length(CurrentClusterMembers[CurrentClusterMembers %in% CurrentUsedRegulators]) != 0){
        # keeping track of the removed Cluster Members:
        CurrentClusterMembers = setdiff(CurrentClusterMembers,CurrentUsedRegulators)  #problem here if the current cluster is empty
        nrCurrentClusterMembers = length(CurrentClusterMembers)
        
        if (nrCurrentClusterMembers > 0){
          names = Data_rownames %in% CurrentClusterMembers
          if (length(which(names==TRUE))>1){
            y = apply((Data[names,]),2,mean) # only removing the used regulators from
          } else {
            y = Data[names,]
          }
          
          fit = cv.glmnet(t(X), y,alpha = alpha, pmax = pmax)
          nonZeroLambdas <- fit$lambda[which(fit$nzero>0)]
          nonZeroCVMs <- fit$cvm[which(fit$nzero>0)]
          
          if(length(which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE)))==0){
            
            #for now: just print a warning, *although* this error WILL cause Amaretto to crash in a few steps.
            warnMessage <- paste0("\nOn cluster ",i," there were no cv.glm results that gave non-zero coefficients during the Autoregulation step.")
            warning(warnMessage);
            
          }
          
          bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE))]
          new_b_o = coef(fit,s = bestNonZeroLambda)
          #what was used up until 09/09/2014 instead of bestNonZeroLambda:
          #new_b_o = coef(fit,s = fit$lambda.1se)
          new_b_opt <- c(new_b_o[2:length(b_o)])
          
          CurrentUsedRegulators = RegulatorData_rownames[which(new_b_opt != 0)]
          nrIterations = nrIterations + 1
          b_opt = new_b_opt
        } else{
          # no more cluster members left, the cluster is empty.
          b_opt = rep(0,length(RegulatorData_rownames))
        }
      }
      Report <- c(length(CurrentClusterPositions),length(CurrentClusterMembers),nrIterations)
      
      #Report(1)=length(CurrentClusterPositions) #original clusters members
      #Report(2)=length(CurrentClusterMembers) #eventual nr of cluster members after removing members that were selected as regulator
      #Report(3)=nrIterations
      #AutoRegulationReport[i,]=t(Report)
      #AutoRegulationReport=Report
    }
    
    # need to do this after the autoregulation, otherwise autoregulation
    # can still add positive microRNA regulators
    if (sum(RegulatorSign[which(RegulatorSign != 0)]) > 0){ # there are limitations on the sign of the regulators
      RegulatorCheck = RegulatorSign * t(b_opt)
      WrongRegulators = which(RegulatorCheck < 0)
      if (length(WrongRegulators)  == 0){# just remove the wrong regulators
        b_opt[WrongRegulators] = 0
      }
    }
    
    if (AutoRegulation >= 1){
      #Beta[i,] = b_opt
    } else {
      BetaSpecial[i] = b_opt
      RegulatorPositions[i] = (RegulatorData_rownames %in% setdiff(RegulatorData_rownames,Data_rownames[CurrentClusterPositions])) # keeping track of the regulators' positions
    }
    
    #y_all[i,] = y
    list(b_opt,y,Report, y_std)
  }
  
  #ptm1<- proc.time() - ptm1
  #printf("Elapsed time is %f seconds\n",ptm1[3])
  if (AutoRegulation == 0){
    for (i in 1:NrClusters){
      Beta[i,RegulatorPositions[i]] = BetaSpecial[i]
    }
  }
  
  tmpPos=NrClusters+1
  
  Beta <- do.call(cbind, BetaY_all[1,2:tmpPos])
  Beta = t(Beta);
  colnames(Beta)=RegulatorData_rownames
  rownames(Beta)=gsub('result.','Module_',rownames(Beta))
  
  y_all<-do.call(cbind, BetaY_all[2,2:tmpPos])
  y_all = t(y_all);
  rownames(y_all)=gsub('result.','Module_',rownames(y_all))
  
  AutoRegulationReport<-do.call(cbind, BetaY_all[3,2:tmpPos])
  AutoRegulationReport = t(AutoRegulationReport)
  rownames(AutoRegulationReport)=gsub('result.','Module_',rownames(AutoRegulationReport))
  
  y_std<-cbind(BetaY_all[4,])
  
  # calculating the error
  error = y_all - (Beta %*% RegulatorData)
  result <- list(y_std = y_std, Beta = Beta,error = error,AutoRegulationReport = AutoRegulationReport)
  return(result)
}

AMARETTO_EvaluateTestSet_lincs <- function(AMARETTOresults,MA_Data_TestSet,RegulatorData_TestSet) {
  nrSamples = ncol(MA_Data_TestSet)
  RegulatorNames=rownames(RegulatorData_TestSet)
  
  #Iterating over the Modules
  stats = mat.or.vec(AMARETTOresults$NrModules,14)
  Rsquare = mat.or.vec(AMARETTOresults$NrModules,1)
  RsquareAjusted = mat.or.vec(AMARETTOresults$NrModules,1)
  modules <- list()
  
  for (i in 1:AMARETTOresults$NrModules){
    #check regulator presence
    currentRegulators = RegulatorNames[which(AMARETTOresults$RegulatoryPrograms[i,] != 0)]
    nrPresentRegulators = sum((rownames(RegulatorData_TestSet)  %in% currentRegulators))
    currentPresentRegulators = (currentRegulators %in% rownames(RegulatorData_TestSet))
    stats[i,1] = nrPresentRegulators
    stats[i,2] = length(currentRegulators)
    
    #checking the presence of the clusters
    currentClusterGenes = AMARETTOresults$AllGenes[which(AMARETTOresults$ModuleMembership[,1] == i)]
    nrPresentClusterGenes = sum((rownames(MA_Data_TestSet) %in% currentClusterGenes))
    stats[i,3] = nrPresentClusterGenes
    stats[i,4] = length(currentClusterGenes)
    stats[i,5] = stats[i,3] / stats[i,4] * 100
    
    #predict cluster expression in test set, always calculate but report
    #the totel percentage weight that is represented
    currentWeights = AMARETTOresults$RegulatoryPrograms[i,which(AMARETTOresults$RegulatoryPrograms[i,] != 0)]
    totalWeights = sum(abs(currentWeights))
    presentWeights = currentWeights[currentPresentRegulators]
    presentRegulators = currentRegulators[currentPresentRegulators]
    totalWeightsPercentage = sum(abs(presentWeights)) / totalWeights * 100
    
    modules[[i]] = currentClusterGenes[currentClusterGenes %in% rownames(MA_Data_TestSet)]
    
    # drop=FALSE, this solves the problem when you have only one regulator, so the previous version is not needed.
    if (nrPresentRegulators > 0) {
      predictions = (t(RegulatorData_TestSet[presentRegulators,,drop=FALSE])) %*% (presentWeights) # need to make sure that the first argument remains a matrix.
      predictions = data.matrix(predictions)
      if (length(modules[[i]]) !=0) {
        if (length(currentClusterGenes)>1){
          outcome = apply(MA_Data_TestSet[currentClusterGenes,],2,mean)
          cx <- t(scale(t(MA_Data_TestSet[currentClusterGenes,])))
          module_PCA = svd(cx)
          outcome_PC <- module_PCA$v[,1]
        } else {
          outcome = MA_Data_TestSet[currentClusterGenes,]
        }
        if(cor(predictions,outcome_PC)<0)
        {
          predictions <- -predictions
        }
        predictions<-predictions/sqrt(sum(predictions^2))
        
        residuals = predictions-outcome + mean(outcome)
        residuals_PC = predictions-outcome_PC + mean(outcome_PC)
        
        module_data<-MA_Data_TestSet[currentClusterGenes,]
        inmodule_corr<-abs(cor(t(module_data),predictions))
        
        # print R2 w.r.t PC
        SStot = sum((outcome_PC-mean(outcome_PC))^2)
        SSres = sum(residuals_PC^2)
        RsquarePC = 1 - (SSres / SStot)
        RsquareAjustedPC = RsquarePC - (nrPresentRegulators/(nrSamples - 1 - nrPresentRegulators))*(1-RsquarePC)
        
        
        # using explained variance as metric, since mean square error is not
        # enough, no baseline interpretation possible
        # SSreg=sumsqr(predictions-mean(outcome))
        
        SStot = sum((outcome-mean(outcome))^2)
        #SSres = sum((predictions-outcome)^2)
        SSres = sum(residuals^2)
        
        Rsquare[i] = 1 - (SSres / SStot)
        #Rsquare[i] = SSres
        RsquareAjusted[i] = Rsquare[i] - (nrPresentRegulators/(nrSamples - 1 - nrPresentRegulators))*(1-Rsquare[i])
        MSE = (1/nrSamples) * sum((residuals)^2)
        
        
        stats[i,6] = totalWeightsPercentage
        stats[i,7] = Rsquare[i]
        stats[i,8] = RsquareAjusted[i]
        stats[i,9] = MSE
        stats[i,10] = median(inmodule_corr)
        stats[i,11] = mad(inmodule_corr)
        stats[i,12] = abs(cor(predictions, outcome_PC))
        stats[i,13] = RsquarePC
        stats[i,14] = RsquareAjustedPC
      } else {
        stats[i,6] = totalWeightsPercentage
        stats[i,7] = 0
        stats[i,8] = 0
        stats[i,9] = 0
        stats[i,10] = 0
        stats[i,11] = 0
        stats[i,12] = 0
        stats[i,13] = 0
        stats[i,14] = 0
      }
    } else {
      stats[i,6] = 0
      stats[i,7] = 0
      stats[i,8] = 0
      stats[i,9] = 0
      stats[i,10] = 0
      stats[i,11] = 0
      stats[i,12] = 0
      stats[i,13] = 0
      stats[i,14] = 0
    }
  }
  #stats = list(stats,CellArrayOfNumToCellArrayofString(num2cell(1:AMARETTOresults.N)),{'nrPresReg' 'nrTotalReg' 'nrPresGen' 'nrTotGen' 'percPresGen' 'percWeightPresent' 'Rsquare' 'RsquareAdjusted'})
  
  dimnames(stats) <- list(rownames(stats, do.NULL = FALSE, prefix = "Module_"),
                          c("nrPresReg" ,"nrTotalReg", "nrPresGen", "nrTotGen", "percPresGen",
                            "percWeightPresent", "Rsquare", "RsquareAdjusted","MSE", "inModuleCorMedian", "inModuleCorMAD", "CorPC", "R2PC", "R2AdjPC"))
  #res <- list(stats = stats,modules = modules)
  
  return(stats)
}




AMARETTO_plot_lincs <- function(AMARETTOinit,AMARETTOresults,ModuleNr) {
  # getting the data
  if (ModuleNr>AMARETTOresults$NrModules){
    stop('\tCannot plot Module',ModuleNr,'since the total number of modules is',AMARETTOresults$N,'.\n')
  }
  ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
  currentRegulators = AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
  RegulatorData=AMARETTOinit$RegulatorData[currentRegulators,]
  ModuleGenes=rownames(ModuleData)
  cat('Module',ModuleNr,'has',length(rownames(ModuleData)),'genes and',length(currentRegulators),'regulators for',length(colnames(ModuleData)),'samples.\n')
  
  # Clustering the module itself
  SampleClustering=hclust(dist(t(ModuleData)), method = "complete", members = NULL)    
  GeneClustering=hclust(dist(ModuleData), method = "complete", members = NULL)
  ClustRegulatorData <- RegulatorData[,SampleClustering$order]
  ClustModuleData <- ModuleData[GeneClustering$order,SampleClustering$order]
  ClustCombinedData <- rbind(ClustModuleData,ClustRegulatorData)
  
  # plotting
  heatmap(ClustCombinedData, name = "Gene expression", column_title = paste('Module',ModuleNr), cluster_rows=FALSE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,row_names_gp=gpar(col=c(rep("white",nrow(ModuleData)),rep("black",nrow(RegulatorData))),fontsize=10),
                     column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Module Genes",nrow(ModuleData)),rep("Regulators",nrow(RegulatorData))),gap = unit(5, "mm"),
                     #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
                     col=colorRamp2(c(-max(abs(ClustCombinedData)), 0, max(abs(ClustCombinedData))), c("green", "black", "red")),heatmap_legend_param = list(color_bar = "continuous"))

  
  
  
}
