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

LINKER_compute_graph_all_LM<-function(Data, lincs_idx, pc_idx)
{
  GEA_per_linc<-list()
  i<-1
  
  X<-Data[lincs_idx,]
  
  driverMat<-matrix(data = 0, nrow = length(pc_idx), ncol = length(lincs_idx))
  
  #compute the LM
  for(idx_gene in 1:length(pc_idx))
  {
    y<-Data[pc_idx[idx_gene],]
    for(i in 1:length(lincs_idx))
    {
      x<-t(X)[,i]
      fit = lm(y~x)
      s<-summary(fit)
      driverMat[idx_gene,i]<-s$coefficients[2,"Pr(>|t|)"]<0.05/(length(pc_idx)*length(lincs_idx))
    }
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

NET_networks_creation<-function(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
{
  all_g<-list()
  all_g$VBSR<-NET_compute_graph_all_VBSR(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
  all_g$LASSOmin<-NET_compute_graph_all_LASSOmin(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
  all_g$LASSO1se<-NET_compute_graph_all_LASSO1se(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
  all_g$LM<-NET_compute_graph_all_LM(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
  return(all_g)
}