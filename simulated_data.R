generate_sim_data<-function(regulatorData){
  
  NrSamples<-ncol(regulatorData)
  NrLincs<-nrow(regulatorData)
  NrModules<-40
  
  maxNrGenes<-10000
  
  ExpectedNrGenesPerModule<-100
  ExpectedNrLincsPerModule<-5
  
  NOISE_SD<-3
  
  Data<-matrix(data=NA, nrow = maxNrGenes, ncol = NrSamples )
  #samples<-rnorm(NrLincs*NrSamples, mean= 0, sd = 1)
  #lincMatrix<-matrix(samples, nrow=NrLincs, ncol = NrSamples)
  
  modulesLincs<-matrix(data=NA,nrow=NrLincs, ncol=NrModules)
  
  ClusterFirstRow<-1
  trueClustering<-numeric(length = maxNrGenes)
  for(i in 1:NrModules){
    moduleLincs<-numeric(NrLincs)
    while(sum(moduleLincs)==0){
      moduleLincs<-rbinom(NrLincs, 1, ExpectedNrLincsPerModule/NrLincs)
    }
    regulatorLincs<-which(moduleLincs==1)
    betas<-rnorm(length(regulatorLincs), 0,1)
    moduleLincs[regulatorLincs]<-betas
    modulesLincs[,i]<-moduleLincs
    module_mean<-moduleLincs%*%regulatorData
    module_mean<-scale(as.vector(module_mean))
    
    NrGenesInModule<-rbinom(1, maxNrGenes, ExpectedNrGenesPerModule/maxNrGenes)
    
    gene_noise<-rnorm(NrSamples*NrGenesInModule,mean=0, sd = NOISE_SD)
    noiseMatrix<-matrix(gene_noise, nrow=NrSamples, ncol = NrGenesInModule)
    geneModule<-as.vector(module_mean)+(noiseMatrix)
    #geneModule<-(sample(c(-1,1), size=NrGenesInModule, replace=TRUE, prob=c(0.5,0.5)))*geneModule
    geneModule<-t(geneModule)
    
    if(is.na(geneModule)){
      foo<-1
    }
    trueClustering[ClusterFirstRow:(ClusterFirstRow+NrGenesInModule-1)]<-i
    Data[ClusterFirstRow:(ClusterFirstRow+NrGenesInModule-1),]<-geneModule
    ClusterFirstRow<-ClusterFirstRow+NrGenesInModule
    
  }
  
  trueClustering<-trueClustering[1:(ClusterFirstRow-1)]
  names(trueClustering)<-sapply(1:(ClusterFirstRow-1), toString)
  
  Data[(ClusterFirstRow):(ClusterFirstRow+NrLincs-1),]<-regulatorData
  
  Data<-Data[1:(ClusterFirstRow+NrLincs-1),]
  rownames(Data)<-c(sapply(1:(ClusterFirstRow-1), toString), rownames(regulatorData))
  

  return(list(Data, 1:(ClusterFirstRow-1), (ClusterFirstRow):(ClusterFirstRow+NrLincs-1), modulesLincs, trueClustering))
  
  
}