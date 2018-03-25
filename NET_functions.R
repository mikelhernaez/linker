NET_compute_graph_all_VBSR<-function(Data, lincs_idx, pc_idx)
{
  
  X<-Data[lincs_idx,]
  
  driverMat<-matrix(data = NA, nrow = length(pc_idx), ncol = length(lincs_idx))
  
  #compute the VBSR
  driverMat<-foreach(idx_gene=1:length(pc_idx), .combine = rbind, .packages = "vbsr")%dopar%
  {
    y<-Data[pc_idx[idx_gene],]
    res<-vbsr(y,t(X),n_orderings = 15,family='normal')
    betas<-res$beta
    betas[res$pval > 0.05/(length(pc_idx)*length(lincs_idx))]<-0
    betas
  }
  rownames(driverMat)<-rownames(Data)[pc_idx]
  colnames(driverMat)<-rownames(Data)[lincs_idx]
  
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_lncRNAs<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_lncRNAs]
  
  g<-graph_from_incidence_matrix(driverMat)
  
  return(g)
  
}

NET_compute_graph_all_LASSOmin<-function(Data, lincs_idx, pc_idx, alpha=1-1e-06)
{
  
  X<-Data[lincs_idx,]
  
  driverMat<-matrix(data = NA, nrow = length(pc_idx), ncol = length(lincs_idx))
  
  #compute the LASSOmin
  driverMat<-foreach(idx_gene=1:length(pc_idx), .combine = rbind, .packages="glmnet")%dopar%
  {
    y<-Data[pc_idx[idx_gene],]
    fit = cv.glmnet(t(X), y, alpha = alpha)
    
    b_o = coef(fit,s = fit$lambda.min)
    b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
    b_opt
  }

  rownames(driverMat)<-rownames(Data)[pc_idx]
  colnames(driverMat)<-rownames(Data)[lincs_idx]
  
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_lncRNAs<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_lncRNAs]
  
  g<-graph_from_incidence_matrix(driverMat)
  
  return(g)
  
}

NET_compute_graph_all_LASSO1se<-function(Data, lincs_idx, pc_idx, alpha=1-1e-06)
{
  
  X<-Data[lincs_idx,]
  
  driverMat<-matrix(data = NA, nrow = length(pc_idx), ncol = length(lincs_idx))
  
  #compute the LASSO1se
  driverMat<-foreach(idx_gene=1:length(pc_idx), .combine = rbind, .packages="glmnet")%dopar%
  {
    y<-Data[pc_idx[idx_gene],]
    fit = cv.glmnet(t(X), y, alpha = alpha)
    
    b_o = coef(fit,s = fit$lambda.1se)
    b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
    b_opt
  }
  
  rownames(driverMat)<-rownames(Data)[pc_idx]
  colnames(driverMat)<-rownames(Data)[lincs_idx]
  
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_lncRNAs<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_lncRNAs]
  
  g<-graph_from_incidence_matrix(driverMat)
  
  return(g)
  
}

NET_compute_graph_all_LM<-function(Data, lincs_idx, pc_idx)
{
  GEA_per_linc<-list()
  i<-1
  
  X<-t(Data[lincs_idx,])
  
  
  #driverMat<-matrix(data = 0, nrow = length(pc_idx), ncol = length(lincs_idx))
  
  #pc_idx=1:100
  #lincs_idx=1:10
  #compute the LM
  NrTotalEdges<-length(pc_idx)*length(lincs_idx)
  Pthre<-0.05/(NrTotalEdges)
  driverMat<-foreach(idx_gene=1:length(pc_idx), .combine=rbind)%dopar%
  #for(idx_gene in 1:length(pc_idx))
    
  {
    y<-Data[pc_idx[idx_gene],]
    NET_compute_LM_from_gene(y,X,Pthre)
  }

  protein_coding_genes<-rownames(Data)[pc_idx]
  lncRNAs<-rownames(Data)[lincs_idx]
  
  rownames(driverMat)<-protein_coding_genes
  colnames(driverMat)<-lncRNAs
  
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_lncRNAs<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_lncRNAs]
  
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

NET_compute_linc_enrichment_from_graph<-function(g, Gene_set_Collections,FDR=0.05, BC=1, NrCores=1)
{
  
  linc_neighbors<-sapply(V(g)[V(g)$type==1], function(x) neighbors(g,x))
  
  registerDoParallel(NrCores)
  
  GEA<-mclapply(linc_neighbors,function(x) gene_enrichment_per_module(names(x), Gene_set_Collections, FDR, BC))
  
  path_lincs<-unlist(lapply(GEA,function(x) unlist(x)))
  unique_paths_linc<-unique(names(path_lincs))
  GEA_per_linc<-path_lincs[unique_paths_linc]
  
  return(GEA_per_linc)
  
}

NET_compute_graph_enrichment_geneSets<-function(pathway_genes,g,FDR=0.05,BC=1, NrCores=1)
{
  GEA<-list()
  Num_lincs<-sum(V(g)$type==1)
  #BC<-Num_lincs*BC
  # BIOCARTA
  Gene_set_Collections<-pathway_genes[c(1)]
  GEA$BIOCARTA<-NET_compute_linc_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC, NrCores=NrCores)
  # KEGG
  Gene_set_Collections<-pathway_genes[c(2)]
  GEA$KEGG<-NET_compute_linc_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC,NrCores=NrCores)
  # REACTOME
  Gene_set_Collections<-pathway_genes[c(3)]
  GEA$REACTOME<-NET_compute_linc_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC,NrCores=NrCores)
  # GENESIGDB
  Gene_set_Collections<-pathway_genes[c(4)]
  GEA$GENESIGDB<-NET_compute_linc_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC,NrCores=NrCores)
  # ALL
  #Gene_set_Collections<-pathway_genes[c(3,4,5,12)]
  #GEA$ALL<-LINKER_compute_linc_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC)
  
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

NET_networks_creation<-function(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx,pathway_genes)
{
  all_g<-list()
  all_g$VBSR<-NET_compute_graph_all_VBSR(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
  print("VBSR done!")
  all_g$LASSOmin<-NET_compute_graph_all_LASSOmin(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
  print("LASSOmin done!")
  all_g$LASSO1se<-NET_compute_graph_all_LASSO1se(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
  print("LASSO1se done!")
  all_g$LM<-NET_compute_graph_all_LM(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
  print("LM done!")
  
  all_GEAs<-NET_compute_graph_enrichment_geneSets_graph_list(pathway_genes,all_g,FDR=0.05,BC=1)
    
  return(list(graphs=all_g,GEAs=all_GEAs))
}