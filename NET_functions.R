NET_run<-function(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx, Gene_set_Collections,
                  graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                  FDR=0.05,
                  NrCores=30)
{
  graphs<-list()
  for(j in 1:length(graph_mode)){
    graphs[[ graph_mode[j] ]] <- switch( graph_mode[j],
                                         "VBSR" = NET_compute_graph_all_VBSR(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx),
                                         "LASSOmin" = NET_compute_graph_all_LASSOmin(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx),
                                         "LASSO1se" = NET_compute_graph_all_LASSO1se(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx),
                                         "LM" = NET_compute_graph_all_LM(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx)
    )
    print(paste0("Graphs for (" ,graph_mode[j], ") computed!"))    
  }
  
  all_GEAs<-NET_compute_graph_enrichment_geneSets_graph_list(pathway_genes,graphs,FDR=FDR,BC=1, NrCores = NrCores)
  
  return(list(graphs=graphs,GEAs=all_GEAs))
}

NET_compute_graph_all_VBSR<-function(Data, regulators_idx, target_idx)
{
  
  X<-Data[regulators_idx,]
  
  driverMat<-matrix(data = NA, nrow = length(target_idx), ncol = length(regulators_idx))
  
  #compute the VBSR
  driverMat<-foreach(idx_gene=1:length(target_idx), .combine = rbind, .packages = "vbsr")%dopar%
  {
    y<-Data[target_idx[idx_gene],]
    res<-vbsr(y,t(X),n_orderings = 15,family='normal')
    betas<-res$beta
    betas[res$pval > 0.05/(length(target_idx)*length(regulators_idx))]<-0
    betas
  }
  rownames(driverMat)<-rownames(Data)[target_idx]
  colnames(driverMat)<-rownames(Data)[regulators_idx]
  
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]
  
  g<-graph_from_incidence_matrix(driverMat)
  
  return(g)
  
}

NET_compute_graph_all_LASSOmin<-function(Data, regulators_idx, target_idx, alpha=1-1e-06)
{
  
  X<-Data[regulators_idx,]
  
  driverMat<-matrix(data = NA, nrow = length(target_idx), ncol = length(regulators_idx))
  
  #compute the LASSOmin
  driverMat<-foreach(idx_gene=1:length(target_idx), .combine = rbind, .packages="glmnet")%dopar%
  {
    y<-Data[target_idx[idx_gene],]
    fit = cv.glmnet(t(X), y, alpha = alpha)
    
    b_o = coef(fit,s = fit$lambda.min)
    b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
    b_opt
  }

  rownames(driverMat)<-rownames(Data)[target_idx]
  colnames(driverMat)<-rownames(Data)[regulators_idx]
  
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]
  
  g<-graph_from_incidence_matrix(driverMat)
  
  return(g)
  
}

NET_compute_graph_all_LASSO1se<-function(Data, regulators_idx, target_idx, alpha=1-1e-06)
{
  
  X<-Data[regulators_idx,]
  
  driverMat<-matrix(data = NA, nrow = length(target_idx), ncol = length(regulators_idx))
  
  #compute the LASSO1se
  driverMat<-foreach(idx_gene=1:length(target_idx), .combine = rbind, .packages="glmnet")%dopar%
  {
    y<-Data[target_idx[idx_gene],]
    fit = cv.glmnet(t(X), y, alpha = alpha)
    
    b_o = coef(fit,s = fit$lambda.1se)
    b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
    b_opt
  }
  
  rownames(driverMat)<-rownames(Data)[target_idx]
  colnames(driverMat)<-rownames(Data)[regulators_idx]
  
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]
  
  g<-graph_from_incidence_matrix(driverMat)
  
  return(g)
  
}

NET_compute_graph_all_LM<-function(Data, regulators_idx, target_idx)
{
  GEA_per_regulator<-list()
  i<-1
  
  X<-t(Data[regulators_idx,])
  
  
  #driverMat<-matrix(data = 0, nrow = length(target_idx), ncol = length(regulators_idx))
  
  #target_idx=1:100
  #regulators_idx=1:10
  #compute the LM
  NrTotalEdges<-length(target_idx)*length(regulators_idx)
  Pthre<-0.05/(NrTotalEdges)
  driverMat<-foreach(idx_gene=1:length(target_idx), .combine=rbind)%dopar%
  #for(idx_gene in 1:length(target_idx))
    
  {
    y<-Data[target_idx[idx_gene],]
    #NET_compute_LM_from_gene(y,X,Pthre)
    driverVec<-numeric(length=ncol(X))
    for(i in 1:ncol(X))
    {
      x<-X[,i]
      fit = lm(y~x)
      s<-summary(fit)
      driverVec[i]<-(s$coefficients[2,"Pr(>|t|)"] < Pthre)
    }
    driverVec
  }

  target_genes<-rownames(Data)[target_idx]
  regulators<-rownames(Data)[regulators_idx]
  
  rownames(driverMat)<-target_genes
  colnames(driverMat)<-regulatory_genes
  
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]
  
  g<-graph_from_incidence_matrix(driverMat)
  
  return(g)
  
}

NET_compute_LM_from_gene<-function(y,X,Pthre){
  
  driverVec<-numeric(length=ncol(X))
  for(i in 1:ncol(X))
  {
    x<-X[,i]
    fit = lm(y~x)
    s<-summary(fit)
    driverVec[i]<-(s$coefficients[2,"Pr(>|t|)"] < Pthre)
  }
  return(driverVec)
}

##################### Enrichment Functions ##############################
NET_compute_regulator_enrichment_from_graph<-function(g, Gene_set_Collections,FDR=0.05, BC=1, NrCores=1)
{
  
  regulator_neighbors<-sapply(V(g)[V(g)$type==1], function(x) neighbors(g,x))
  
  registerDoParallel(NrCores)
  
  GEA<-mclapply(regulator_neighbors,function(x) NET_module_gene_enrichment(names(x), Gene_set_Collections, FDR, BC))
  
  path_regulators<-unlist(lapply(GEA,function(x) unlist(x)))
  unique_paths_regulator<-unique(names(path_regulators))
  GEA_per_regulator<-path_regulators[unique_paths_regulator]
  
  return(GEA_per_regulator)
  
}

NET_compute_graph_enrichment_geneSets<-function(pathway_genes,g,FDR=0.05,BC=1, NrCores=1)
{
  GEA<-list()
  Num_regulators<-sum(V(g)$type==1)
  #BC<-Num_regulators*BC
  # BIOCARTA
  Gene_set_Collections<-pathway_genes[c(1)]
  GEA$BIOCARTA<-NET_compute_regulator_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC, NrCores=NrCores)
  # KEGG
  Gene_set_Collections<-pathway_genes[c(2)]
  GEA$KEGG<-NET_compute_regulator_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC,NrCores=NrCores)
  # REACTOME
  Gene_set_Collections<-pathway_genes[c(3)]
  GEA$REACTOME<-NET_compute_regulator_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC,NrCores=NrCores)
  # GENESIGDB
  Gene_set_Collections<-pathway_genes[c(4)]
  GEA$GENESIGDB<-NET_compute_regulator_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC,NrCores=NrCores)
  # ALL
  #Gene_set_Collections<-pathway_genes[c(3,4,5,12)]
  #GEA$ALL<-LINKER_compute_regulator_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC)
  
  return(GEA)
}

NET_compute_graph_enrichment_geneSets_graph_list<-function(pathway_genes,g,FDR=0.05,BC=1, NrCores=1)
{
  
  GEA<-list()
  
  for(i in 1:length(g))
  {
    net_mode<-names(g)[i]

    GEA[[ net_mode ]]<-NET_compute_graph_enrichment_geneSets(pathway_genes,g[[net_mode]], BC=BC, NrCores=NrCores)
    print(paste0("Mode ",net_mode," computed!"))
  }
  
  return(GEA)
}

######################### Main Enrichment Test Function ########################
NET_module_gene_enrichment <-function(Module, pathway_genes, FDR=0.05, BC=1, total_genes=10000)
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

