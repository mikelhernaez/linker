generate_enrichment_file_per_module<-function(Module_number, results,  Gene_set_Collection, file_name, FDR)
{
  #Module_protein_coding_genes<-names(which(results[[3]]$ModuleMembership[,]==Module_number))
  Module_protein_coding_genes<-results$AllGenes[which(results$ModuleMembership[,]==Module_number)]
  Module_protein_coding_gene_list<-sapply(Module_protein_coding_genes, function(x) strsplit(x, "\\|"))
  Module_protein_coding_genes<-sapply(Module_protein_coding_gene_list, function(x) x[[6]])
  Module_protein_coding_genes<-unname(Module_protein_coding_genes)
  
  Modules_lincRNAs<-names(which(results$RegulatoryPrograms[Module_number,]!=0))
  Modules_lincRNAs_list<-sapply(Modules_lincRNAs, function(x) strsplit(x, "\\|"))
  Modules_lincRNAs<-sapply(Modules_lincRNAs_list, function(x) x[[6]])
  Modules_lincRNAs<-unname(Modules_lincRNAs)
  
  
  GEA<-gene_enrichment_per_module(Module_protein_coding_genes, Gene_set_Collection, FDR)
  
  if(sum(sapply(GEA, length))==0){
    return()
  }
  
  
  
  # write(paste(">Module",Module_number, collapse = " "), file_name)
  # write.table(results[[1]][Module_number,14], file_name, append = TRUE)
  # 
  # write(">lincRNAs", file_name, sep=" ", append = TRUE)
  # write.table(Modules_lincRNAs, file_name, sep="\t", append = TRUE)
  # write(">protein_coding_genes", file_name, sep=" ", append = TRUE)
  # write.table(Module_protein_coding_genes, file_name,sep="\t" ,append=TRUE)
  # 
  # for(i in 1:length(GEA)){
  #   write(paste(">",i, collapse = ""), file_name, append = TRUE)
  #   df_gea<-as.data.frame(GEA[[i]])
  #   write.table(df_gea,file_name, sep="\t",append=TRUE)
  # }
}





gene_enrichment_per_module <-function(Module, pathway_genes, FDR, BC=1)
{

  total_genes<-10000
  
  i<-1
  under_enrichment_pvalues<-numeric()
  over_enrichment_pvalues<-numeric()
  for(collection_idx in 1:length(pathway_genes)){
  
    for(pathway_idx in 1:length(pathway_genes[[collection_idx]])){
    
      white_balls<-length(pathway_genes[[collection_idx]][[pathway_idx]])
      black_balls<-total_genes - white_balls
      
      drawn<-length(Module)
      drawn_whites<-length(which(pathway_genes[[collection_idx]][[pathway_idx]] %in% Module))
      
      under_enrichment_pvalues[i]<- phyper(drawn_whites, white_balls, black_balls, drawn, lower.tail = TRUE, log.p = FALSE)
      over_enrichment_pvalues[i]<- phyper(drawn_whites-1, white_balls, black_balls, drawn, lower.tail = FALSE, log.p = FALSE)
      i<-i+1
    }
  }

  if(BC<0){
    under_adjp<-under_adjp*(-BC)
    over_adjp<-over_adjp*(-BC)
  }
  
  else{
    under_adjp<-p.adjust(under_enrichment_pvalues,'fdr')
    over_adjp<-p.adjust(over_enrichment_pvalues,'fdr')
    under_adjp<-under_adjp*BC
    over_adjp<-over_adjp*BC
  }
      
      

  i<-1
  GEA_over<-list()
  GEA_under<-list()
  for(collection_idx in 1:length(pathway_genes))  
  {
    j<-1
    jj<-1
    GEA_under[[collection_idx]]<-numeric()
    GEA_over[[collection_idx]]<-numeric()
    gensetNames_over<-character()
    gensetNames_under<-character()
    for(pathway_idx in 1:length(pathway_genes[[collection_idx]]))
    {
      if(over_adjp[i]<FDR){
        GEA_over[[collection_idx]][j]<-over_adjp[i]
        gensetNames_over[j]<-names(pathway_genes[[collection_idx]])[pathway_idx]
        j<-j+1
      }
      
      if(under_adjp[i]<FDR){
        GEA_under[[collection_idx]][jj]<-over_adjp[i]
        gensetNames_under[jj]<-names(pathway_genes[[collection_idx]])[pathway_idx]
        jj<-jj+1
      }
      i<-i+1
    }
    names(GEA_over[[collection_idx]])<-gensetNames_over
    names(GEA_under[[collection_idx]])<-gensetNames_under
  }
  #return(list(GEA_under, GEA_over))
  return(GEA_over)
}

filter_enriched_modules<-function(Gene_set_Collections,results,FDR=0.05){
  
  NrSims<-length(results)
  

  enriched_modules<-list()
  enriched_idx<-1
  for(idx_sim in 1:NrSims){
    NrBootstraps<-length(results[[idx_sim]]$bootstrapResult)
    for(idx_bootstrap in 1:NrBootstraps){
      
      NrModules<-results[[idx_sim]]$bootstrapResult[[idx_bootstrap]]$NrModules
      boot_results<-results[[idx_sim]]$bootstrapResult[[idx_bootstrap]]
      for(Module_number in 1: NrModules){
        
        #file_name<-paste(c("test1",module_idx), collapse = '_')
        Module_protein_coding_genes_full_name<-boot_results$AllGenes[which(boot_results$ModuleMembership[,]==Module_number)]
        Module_protein_coding_gene_list<-sapply(Module_protein_coding_genes_full_name, function(x) strsplit(x, "\\|"))
        Module_protein_coding_genes<-sapply(Module_protein_coding_gene_list, function(x) x[[6]])
        Module_protein_coding_genes<-unname(Module_protein_coding_genes)
        #Module_protein_coding_genes<-Module_protein_coding_genes_full_name
        
        Modules_lncRNAs_full_name<-names(which(boot_results$RegulatoryPrograms[Module_number,]!=0))
        Modules_lncRNAs_list<-sapply(Modules_lncRNAs_full_name, function(x) strsplit(x, "\\|"))
        Modules_lncRNAs<-sapply(Modules_lncRNAs_list, function(x) x[[6]])
        Modules_lncRNAs<-unname(Modules_lncRNAs)
        #Modules_lncRNAs<-Modules_lncRNAs_full_name
        
        GEA<-gene_enrichment_per_module(Module_protein_coding_genes, Gene_set_Collections, FDR)
        
        if(sum(sapply(GEA, length))==0){
          #next()
        }
        enriched_modules[[enriched_idx]]<-list(
            GEA=GEA,
            protein_coding_genes=Module_protein_coding_genes_full_name, 
             lncRNAs=Modules_lncRNAs_full_name, 
             regulatory_program=boot_results$RegulatoryPrograms[Module_number,],
             training_stats=boot_results$trainingStats[Module_number,],
             test_stats=results[[idx_sim]]$bootstrapTestStats[[idx_bootstrap]][Module_number],
             assigned_genes=which(boot_results$ModuleMembership[,]==Module_number))
        
        enriched_idx<-enriched_idx + 1
      }
    }
  }
  return(enriched_modules)
  
}

LINKER_plot_enrichement_bootstrap<-function(Gene_set_Collections,results,FDR=0.05){
  
  NrSims<-length(results)
  
  
  enriched_modules<-list()
  enriched_idx<-1
  for(idx_sim in 1:NrSims){
    NrBootstraps<-length(results[[idx_sim]]$bootstrapResult)
    GEA_boot<-list()
    for(idx_bootstrap in 1:NrBootstraps){
      
      NrModules<-results[[idx_sim]]$bootstrapResult[[idx_bootstrap]]$NrModules
      boot_results<-results[[idx_sim]]$bootstrapResult[[idx_bootstrap]]
      GEA_boot[[idx_bootstrap]]<-list()
      for(Module_number in 1: NrModules){
        
        #file_name<-paste(c("test1",module_idx), collapse = '_')
        Module_protein_coding_genes_full_name<-boot_results$AllGenes[which(boot_results$ModuleMembership[,]==Module_number)]
        Module_protein_coding_gene_list<-sapply(Module_protein_coding_genes_full_name, function(x) strsplit(x, "\\|"))
        Module_protein_coding_genes<-sapply(Module_protein_coding_gene_list, function(x) x[[6]])
        Module_protein_coding_genes<-unname(Module_protein_coding_genes)
        #Module_protein_coding_genes<-Module_protein_coding_genes_full_name
        
        Modules_lncRNAs_full_name<-names(which(boot_results$RegulatoryPrograms[Module_number,]!=0))
        Modules_lncRNAs_list<-sapply(Modules_lncRNAs_full_name, function(x) strsplit(x, "\\|"))
        Modules_lncRNAs<-sapply(Modules_lncRNAs_list, function(x) x[[6]])
        Modules_lncRNAs<-unname(Modules_lncRNAs)
        #Modules_lncRNAs<-Modules_lncRNAs_full_name
        
        GEA<-gene_enrichment_per_module(Module_protein_coding_genes, Gene_set_Collections, FDR)
        GEA_boot[[idx_bootstrap]][[Module_number]]<-unlist(GEA)
      }
    }
    
  }
  
}

LINKER_compute_enriched_lincs_communities<-function(enriched_modules, enriched_graphs, Data)
{
  
  GEA_g<-list()
  GEA_com<-list()
  GEA_per_linc<-list()
  bp_g<-list()
  i<-1
  
  graph_community<-enriched_graphs[[2]]
  NrCommunities<-length(unique(graph_community$membership))
  for(comm_idx in 1:NrCommunities){
    community_modules<-enriched_modules[which(graph_community$membership==comm_idx)]
    community_pc_genes<-unique(unlist(lapply(community_modules, function(x) x$protein_coding_genes)))
    community_lncRNAs<-unique(unlist(lapply(community_modules, function(x) x$lncRNAs)))
    community_GEA_paths<-unlist(lapply(community_modules, function(x) x$GEA))
    community_GEA_unique_paths<-unique(names( community_GEA_paths ))
    community_GEA<-community_GEA_paths[community_GEA_unique_paths]
    X<-Data[community_lncRNAs,]
    
    driverMat<-matrix(data = NA, nrow = length(community_pc_genes), ncol = length(community_lncRNAs))
    
    #We need to handle the special case where only one lincRNA regulates a module/community
    if(length(community_lncRNAs)<2){

      non_zero_beta<-community_modules[[mod_idx]]$regulatory_program[which(community_modules[[mod_idx]]$regulatory_program != 0)]
      if(length(non_zero_beta) != 1){
        print("error")
      }
      
      driverMat<-matrix(data = non_zero_beta, nrow = length(community_pc_genes), ncol = length(community_lncRNAs))
      
      Module_protein_coding_gene_list<-sapply(community_pc_genes, function(x) strsplit(x, "\\|"))
      Module_protein_coding_genes<-unname(sapply(Module_protein_coding_gene_list, function(x) x[[6]]))
      Module_lncRNAs_list<-sapply(community_lncRNAs, function(x) strsplit(x, "\\|"))
      Module_lncRNAs<-unname(sapply(Module_lncRNAs_list, function(x) x[[6]]))
      
      rownames(driverMat)<-Module_protein_coding_genes
      colnames(driverMat)<-Module_lncRNAs
    }
    else{
      for(idx_gene in 1:length(community_pc_genes))
      {
        y<-Data[community_pc_genes[idx_gene],]
        #print(paste0("(",idx_gene,",",comm_idx,")",":",length(y),";", dim(X)))
        res<-vbsr(y,t(X),n_orderings = 15,family='normal')
        betas<-res$beta
        betas[res$pval > 0.05/(length(community_lncRNAs)*length(community_pc_genes))]<-0
        driverMat[idx_gene,]<-betas
      }
      
      Module_protein_coding_gene_list<-sapply(community_pc_genes, function(x) strsplit(x, "\\|"))
      Module_protein_coding_genes<-unname(sapply(Module_protein_coding_gene_list, function(x) x[[6]]))
      Module_lncRNAs_list<-sapply(community_lncRNAs, function(x) strsplit(x, "\\|"))
      Module_lncRNAs<-unname(sapply(Module_lncRNAs_list, function(x) x[[6]]))
      
      rownames(driverMat)<-Module_protein_coding_genes
      colnames(driverMat)<-Module_lncRNAs
      
      regulated_genes<-which(rowSums(abs(driverMat))!=0)
      regulatory_lncRNAs<-which(colSums(abs(driverMat))!=0)
      
      # We need to treat the special cases independently
      if(length(regulated_genes)<2){
        driverMat<-driverMat[regulated_genes,]
        driverMat<-driverMat[regulatory_lncRNAs]
      }
      else if(length(regulatory_lncRNAs)<2){
        driverMat<-driverMat[,regulatory_lncRNAs]
        driverMat<-driverMat[regulated_genes]
      }
      else{
        driverMat<-driverMat[,regulatory_lncRNAs]
        driverMat<-driverMat[regulated_genes,]
      }
      

      
    }
    
    g<-graph_from_incidence_matrix(driverMat)
    #g_comps<-decompose.graph(g)
    #g_comps_size<-sapply(g_comps, vcount)
    #g<-g_comps[[which(g_comps_size==max(g_comps_size))]]
    
    
    #non_zero_idx<-which(driverMat!=0, arr.ind = TRUE)
    #sparse_driverMat<-sparseMatrix(non_zero_idx[,1], non_zero_idx[,2], x=driverMat[which(driverMat!=0)])
    
    col <- c("steelblue", "orange")
    shape <- c("circle", "square")
    #plot(g, vertex.label=NA, vertex.size=7, layout=layout_as_tree,vertex.color = col[as.numeric(V(g)$type)+1],
    #     vertex.shape = shape[as.numeric(V(g)$type)+1])
    
    GEA_g_comm<-gene_enrichment_per_module(names(V(g)[V(g)$type==0]), Gene_set_Collections, FDR=0.05)
    
    linc_neighbors<-sapply(V(g)[V(g)$type==1], function(x) neighbors(g,x))
    
    GEA_per_linc_comm<-lapply(linc_neighbors,function(x) gene_enrichment_per_module(names(x), Gene_set_Collections, FDR=0.05))
    
    GEA_g[[i]]<-unlist(GEA_g_comm)
    GEA_com[[i]]<-community_GEA
    path_lincs<-unlist(unname(lapply(GEA_per_linc_comm,function(x) unlist(x))))
    unique_paths_linc<-unique(names(path_lincs))
    GEA_per_linc[[i]]<-path_lincs[unique_paths_linc]
    bp_g[[i]]<-g
    i<-i+1
    
    
    
  }
  
  return(list(GEA_module=GEA_com,GEA_graph=GEA_g, GEA_lincs=GEA_per_linc, graph=bp_g))
  
}

LINKER_compute_modules_graph<-function(modules, Data, mode="VBSR",alpha=1-1e-06)
{
  
  bp_g<-list()
  i<-1
    
  for(mod_idx in 1:length(modules))
  {
    pc_genes<-unlist(modules[[mod_idx]]$protein_coding_genes)
    lncRNAs<-unlist(modules[[mod_idx]]$lncRNAs)
    X<-Data[lncRNAs,]
    
    #We need to handle the special case where only one lincRNA regulates a module/community
    if(length(lncRNAs)<2){
        
      non_zero_beta<-modules[[mod_idx]]$regulatory_program[which(modules[[mod_idx]]$regulatory_program != 0)]
      if(length(non_zero_beta) != 1){
        print("error")
      }
        
      driverMat<-matrix(data = non_zero_beta, nrow = length(pc_genes), ncol = length(lncRNAs))
        
      Module_protein_coding_gene_list<-sapply(pc_genes, function(x) strsplit(x, "\\|"))
      Module_protein_coding_genes<-unname(sapply(Module_protein_coding_gene_list, function(x) x[[6]]))
      Module_lncRNAs_list<-sapply(lncRNAs, function(x) strsplit(x, "\\|"))
      Module_lncRNAs<-unname(sapply(Module_lncRNAs_list, function(x) x[[6]]))
        
      rownames(driverMat)<-Module_protein_coding_genes
      colnames(driverMat)<-Module_lncRNAs
    }
    else{
        
      driverMat<-matrix(data = NA, nrow = length(pc_genes), ncol = length(lncRNAs))
      
      for(idx_gene in 1:length(pc_genes))
      {
        y<-Data[pc_genes[idx_gene],]
        
        if(mode=="VBSR")
        {
          res<-vbsr(y,t(X),n_orderings = 15,family='normal')
          betas<-res$beta
          betas[res$pval > 0.05/(length(lncRNAs)*length(pc_genes))]<-0
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
          for(idx_lincs in 1:length(lncRNAs))
          {
            x<-t(X)[,idx_lincs]
            fit = lm(y~x)
            s<-summary(fit)
            driverMat[idx_gene,idx_lincs]<-s$coefficients[2,"Pr(>|t|)"]<0.05/(length(pc_genes)*length(lncRNAs))
          }
        }
        else
        {
          print("MODE NOT RECOGNIZED")
        }
        
      }
    
      Module_protein_coding_gene_list<-sapply(pc_genes, function(x) strsplit(x, "\\|"))
      Module_protein_coding_genes<-unname(sapply(Module_protein_coding_gene_list, function(x) x[[6]]))
      Module_lncRNAs_list<-sapply(lncRNAs, function(x) strsplit(x, "\\|"))
      Module_lncRNAs<-unname(sapply(Module_lncRNAs_list, function(x) x[[6]]))
    
      rownames(driverMat)<-Module_protein_coding_genes
      colnames(driverMat)<-Module_lncRNAs
    
      regulated_genes<-which(rowSums(abs(driverMat))!=0)
      regulatory_lncRNAs<-which(colSums(abs(driverMat))!=0)
        
      # We need to treat the special cases independently
      if(length(regulated_genes)<2){
        driverMat<-driverMat[regulated_genes,]
        driverMat<-driverMat[regulatory_lncRNAs]
      }
      else if(length(regulatory_lncRNAs)<2){
        driverMat<-driverMat[,regulatory_lncRNAs]
        driverMat<-driverMat[regulated_genes]
      }
      else{
        driverMat<-driverMat[,regulatory_lncRNAs]
        driverMat<-driverMat[regulated_genes,]
      }

    }

    bp_g[[i]]<-graph_from_incidence_matrix(driverMat)
    i<-i+1
  }
    
  return(bp_g)
}



LINKER_compute_graph_all_VBSR<-function(Data, lincs_idx, pc_idx)
{
  
  X<-Data[lincs_idx,]
    
  driverMat<-matrix(data = NA, nrow = length(pc_idx), ncol = length(lincs_idx))
  
  #compute the VBSR
  for(idx_gene in 1:length(pc_idx))
  {
    y<-Data[pc_idx[idx_gene],]
    res<-vbsr(y,t(X),n_orderings = 15,family='normal')
    betas<-res$beta
    betas[res$pval > 0.05/(length(pc_idx)*length(lincs_idx))]<-0
    driverMat[idx_gene,]<-betas
  }
  
  
  protein_coding_gene_list<-sapply(rownames(Data)[pc_idx], function(x) strsplit(x, "\\|"))
  protein_coding_genes<-unname(sapply(protein_coding_gene_list, function(x) x[[6]]))
  lncRNAs_list<-sapply(rownames(Data)[lincs_idx], function(x) strsplit(x, "\\|"))
  lncRNAs<-unname(sapply(lncRNAs_list, function(x) x[[6]]))
    
  rownames(driverMat)<-protein_coding_genes
  colnames(driverMat)<-lncRNAs
    
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_lncRNAs<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_lncRNAs]
    
  g<-graph_from_incidence_matrix(driverMat)
  
  return(g)
  
}

LINKER_compute_graph_all_LASSO_min<-function(Data, lincs_idx, pc_idx, alpha=1-1e-06)
{
  
  X<-Data[lincs_idx,]
  
  driverMat<-matrix(data = NA, nrow = length(pc_idx), ncol = length(lincs_idx))

  #compute the VBSR
  for(idx_gene in 1:length(pc_idx))
  {
    y<-Data[pc_idx[idx_gene],]
    fit = cv.glmnet(t(X), y, alpha = alpha)
    
    b_o = coef(fit,s = fit$lambda.min)
    b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
    driverMat[idx_gene,]<-b_opt
  }

  protein_coding_gene_list<-sapply(rownames(Data)[pc_idx], function(x) strsplit(x, "\\|"))
  protein_coding_genes<-unname(sapply(protein_coding_gene_list, function(x) x[[6]]))
  lncRNAs_list<-sapply(rownames(Data)[lincs_idx], function(x) strsplit(x, "\\|"))
  lncRNAs<-unname(sapply(lncRNAs_list, function(x) x[[6]]))
  
  rownames(driverMat)<-protein_coding_genes
  colnames(driverMat)<-lncRNAs
  
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_lncRNAs<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_lncRNAs]
  
  g<-graph_from_incidence_matrix(driverMat)
  
  return(g)
  
}

LINKER_compute_graph_all_LASSO_1se<-function(Data, lincs_idx, pc_idx, alpha=1-1e-06)
{
  
  X<-Data[lincs_idx,]
  
  driverMat<-matrix(data = NA, nrow = length(pc_idx), ncol = length(lincs_idx))
  
  #compute the VBSR
  for(idx_gene in 1:length(pc_idx))
  {
    y<-Data[pc_idx[idx_gene],]
    fit = cv.glmnet(t(X), y, alpha = alpha)
    
    b_o = coef(fit,s = fit$lambda.1se)
    b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
    driverMat[idx_gene,]<-b_opt
  }
  
  protein_coding_gene_list<-sapply(rownames(Data)[pc_idx], function(x) strsplit(x, "\\|"))
  protein_coding_genes<-unname(sapply(protein_coding_gene_list, function(x) x[[6]]))
  lncRNAs_list<-sapply(rownames(Data)[lincs_idx], function(x) strsplit(x, "\\|"))
  lncRNAs<-unname(sapply(lncRNAs_list, function(x) x[[6]]))
  
  rownames(driverMat)<-protein_coding_genes
  colnames(driverMat)<-lncRNAs
  
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
  
  protein_coding_gene_list<-sapply(rownames(Data)[pc_idx], function(x) strsplit(x, "\\|"))
  protein_coding_genes<-unname(sapply(protein_coding_gene_list, function(x) x[[6]]))
  lncRNAs_list<-sapply(rownames(Data)[lincs_idx], function(x) strsplit(x, "\\|"))
  lncRNAs<-unname(sapply(lncRNAs_list, function(x) x[[6]]))
  
  rownames(driverMat)<-protein_coding_genes
  colnames(driverMat)<-lncRNAs
  
  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_lncRNAs<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_lncRNAs]
  
  g<-graph_from_incidence_matrix(driverMat)
  
  return(g)
  
}
LINKER_compute_linc_enrichment_from_graph<-function(g, Gene_set_Collections,FDR=0.05, BC=1)
{
  
  linc_neighbors<-sapply(V(g)[V(g)$type==1], function(x) neighbors(g,x))
  
  GEA<-lapply(linc_neighbors,function(x) gene_enrichment_per_module(names(x), Gene_set_Collections, FDR, BC))
  
  path_lincs<-unlist(lapply(GEA,function(x) unlist(x)))
  unique_paths_linc<-unique(names(path_lincs))
  GEA_per_linc<-path_lincs[unique_paths_linc]
  
  return(GEA_per_linc)
  
}

LINKER_compute_linc_enrichment_from_graph_list<-function(g_list, Gene_set_Collections,FDR=0.05, BC=1)
{
  
  path_lincs<-list()
  for(i in 1:length(g_list))
  {
    g<-g_list[[i]]
    linc_neighbors<-sapply(V(g)[V(g)$type==1], function(x) neighbors(g,x))
    GEA<-lapply(linc_neighbors,function(x) gene_enrichment_per_module(names(x), Gene_set_Collections, FDR, BC))
    path_lincs[[i]]<-unlist(lapply(GEA,function(x) unlist(x)))
  }
  
  path_lincs<-unlist(path_lincs)
  unique_paths_linc<-unique(names(path_lincs))
  GEA_per_linc<-path_lincs[unique_paths_linc]
  
  return(GEA_per_linc)
  
}

LINKER_compute_graph_enrichment_geneSets<-function(pathway_genes,g,FDR=0.05,BC=1)
{
  GEA<-list()
  Num_lincs<-sum(V(g)$type==1)
  BC<-Num_lincs*BC
  # BIOCARTA
  Gene_set_Collections<-pathway_genes[c(3)]
  GEA$BIOCARTA<-LINKER_compute_linc_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC)
  # KEGG
  Gene_set_Collections<-pathway_genes[c(4)]
  GEA$KEGG<-LINKER_compute_linc_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC)
  # REACTOME
  Gene_set_Collections<-pathway_genes[c(5)]
  GEA$REACTOME<-LINKER_compute_linc_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC)
  # GENESIGDB
  Gene_set_Collections<-pathway_genes[c(12)]
  GEA$GENESIGDB<-LINKER_compute_linc_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC)
  # ALL
  #Gene_set_Collections<-pathway_genes[c(3,4,5,12)]
  #GEA$ALL<-LINKER_compute_linc_enrichment_from_graph(g, Gene_set_Collections,FDR=FDR, BC=BC)
  
  return(GEA)
}

LINKER_compute_graph_list_enrichment_geneSets<-function(pathway_genes,g_list,FDR=0.05,BC=20)
{
  GEA<-list()
  
  Num_lincs<-sapply(g_list, function(x) sum(V(g)$type==1))
  BC<-mean(Num_lincs)*BC
  # BIOCARTA
  Gene_set_Collections<-pathway_genes[c(3)]
  GEA$BIOCARTA<-LINKER_compute_linc_enrichment_from_graph_list(g_list, Gene_set_Collections,FDR=FDR, BC=BC)
  # KEGG
  Gene_set_Collections<-pathway_genes[c(4)]
  GEA$KEGG<-LINKER_compute_linc_enrichment_from_graph_list(g_list, Gene_set_Collections,FDR=FDR, BC=BC)
  # REACTOME
  Gene_set_Collections<-pathway_genes[c(5)]
  GEA$REACTOME<-LINKER_compute_linc_enrichment_from_graph_list(g_list, Gene_set_Collections,FDR=FDR, BC=BC)
  # GENESIGDB
  Gene_set_Collections<-pathway_genes[c(12)]
  GEA$GENESIGDB<-LINKER_compute_linc_enrichment_from_graph_list(g_list, Gene_set_Collections,FDR=FDR, BC=BC)
  # ALL
  #Gene_set_Collections<-pathway_genes[c(3,4,5,12)]
  #GEA$ALL<-LINKER_compute_linc_enrichment_from_graph_list(g_list, Gene_set_Collections,FDR=FDR, BC=BC)
  
  return(GEA)
}

LINKER_compute_graph_enrichment_geneSets_graph_list<-function(pathway_genes,graphs,FDR=0.05,BC=1)
{
  graph_names<-names(graphs)
  
  GEA<-list()
  for(i in 1:length(graphs))
  {
    if(class(graphs[[i]])=="list"){
      
      GEA[[ graph_names[i] ]]<-LINKER_compute_graph_list_enrichment_geneSets(pathway_genes,graphs[[i]], BC=20)
    }
    else{
      GEA[[ graph_names[i] ]]<-LINKER_compute_graph_enrichment_geneSets(pathway_genes,graphs[[i]])
    }
  }
  
  return(GEA)
}

LINKER_plot_graphs_topology<-function(graphs)
{
  L<-length(graphs)
  col=c("black","red","blue","green","black","red","blue","green","black","red","blue","green")
  pch=c(0,2,5,6,18,19,20,1,3,4,7,8)
  
  
  par(mfrow=c(3,4))
  for(idx in 1:length(graphs)){
    if(class(graphs[[idx]])=="igraph"){
      t<-table(degree(graphs[[idx]]))
    }
    else{
      t<-table(unlist(lapply(graphs[[idx]], function(x) degree(x))))
    }
    plot(as.numeric(names(t)),as.numeric(t)/sum(as.numeric(t)), col=col[idx], 
             pch = 16,log="xy",xlab='Degree',ylab='P(Degree)', main = names(graphs)[idx],xlim=c(1,5000), ylim=c(1/500000,1))
  }
  
  par(mfrow=c(3,4))
  for(idx in 1:length(graphs)){
    if(class(graphs[[idx]])=="igraph"){
      t<-table(degree(graphs[[idx]], V(graphs[[idx]])[V(graphs[[idx]])$type==0]))
    }
    else{
      t<-table(unlist(unname(lapply(graphs[[idx]], function(x) degree(x, V(x)[V(x)$type==0] )))))
    }
    plot(as.numeric(names(t)),as.numeric(t)/sum(as.numeric(t)), col=col[idx], 
         pch = 16,log="xy",xlab='Degree',ylab='P(Degree)', main = names(graphs)[idx],xlim=c(1,500), ylim=c(1/500000,1))
  }
  
  par(mfrow=c(3,4))
  for(idx in 1:length(graphs)){
    if(class(graphs[[idx]])=="igraph"){
      t<-table(degree(graphs[[idx]], V(graphs[[idx]])[V(graphs[[idx]])$type==1]))
    }
    else{
      t<-table(unlist(unname(lapply(graphs[[idx]], function(x) degree(x, V(x)[V(x)$type==1] )))))
    }
    plot(as.numeric(names(t)),as.numeric(t)/sum(as.numeric(t)), col=col[idx], 
         pch = 16,log="xy",xlab='Degree',ylab='P(Degree)', main = names(graphs)[idx],xlim=c(1,5000), ylim=c(1/500000,1))
  }
  
}
  
LINKER_plot_GEAs<-function(GEAs)
{  
  
  par()
  plot(1, type="n", xlab="-log(p-value)", ylab="cummulative counts of number of enriched gene sets", xlim=c(-10, -1), ylim=c(0, 50), main="BIOCARTA")
  
  h<-hist(log(GEAs$VBSR$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black")
  
  h<-hist(log(GEAs$LASSOmin$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red")
  
  h<-hist(log(GEAs$LASSO1se$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue")
  
  h<-hist(log(GEAs$LM$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green")
  
  h<-hist(log(GEAs$VBSR_VBSR$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black", type = "b")
  
  h<-hist(log(GEAs$VBSR_LASSOmin$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red", type = "b")
  
  h<-hist(log(GEAs$VBSR_LASSO1se$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue", type = "b")
  
  h<-hist(log(GEAs$VBSR_LM$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green", type = "b")
  
  h<-hist(log(GEAs$LASSO_VBSR$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black", type = "o")
  
  h<-hist(log(GEAs$LASSO_LASSOmin$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red", type = "o")
  
  h<-hist(log(GEAs$LASSO_LASSO1se$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue", type = "o")
  
  h<-hist(log(GEAs$LASSO_LM$BIOCARTA), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green", type = "o")
  
  plot(1, type="n", xlab="-log(p-value)", ylab="cummulative counts of number of enriched gene sets", xlim=c(-15, -1), ylim=c(0, 200), main="KEGG")
  
  h<-hist(log(GEAs$VBSR$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black")
  
  h<-hist(log(GEAs$LASSOmin$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red")
  
  h<-hist(log(GEAs$LASSO1se$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue")
  
  h<-hist(log(GEAs$LM$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green")
  
  h<-hist(log(GEAs$VBSR_VBSR$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black", type = "b")
  
  h<-hist(log(GEAs$VBSR_LASSOmin$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red", type = "b")
  
  h<-hist(log(GEAs$VBSR_LASSO1se$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue", type = "b")
  
  h<-hist(log(GEAs$VBSR_LM$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green", type = "b")
  
  h<-hist(log(GEAs$LASSO_VBSR$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black", type = "o")
  
  h<-hist(log(GEAs$LASSO_LASSOmin$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red", type = "o")
  
  h<-hist(log(GEAs$LASSO_LASSO1se$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue", type = "o")
  
  h<-hist(log(GEAs$LASSO_LM$KEGG), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green", type = "o")
  
  
  
  
  plot(1, type="n", xlab="-log(p-value)", ylab="cummulative counts of number of enriched gene sets", xlim=c(-20, -1), ylim=c(0, 200), main="REACTOME")
  
  h<-hist(log(GEAs$VBSR$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black")
  
  h<-hist(log(GEAs$LASSOmin$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red")
  
  h<-hist(log(GEAs$LASSO1se$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue")
  
  h<-hist(log(GEAs$LM$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green")
  
  h<-hist(log(GEAs$VBSR_VBSR$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black", type = "b")
  
  h<-hist(log(GEAs$VBSR_LASSOmin$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red", type = "b")
  
  h<-hist(log(GEAs$VBSR_LASSO1se$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue", type = "b")
  
  h<-hist(log(GEAs$VBSR_LM$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green", type = "b")
  
  h<-hist(log(GEAs$LASSO_VBSR$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black", type = "o")
  
  h<-hist(log(GEAs$LASSO_LASSOmin$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red", type = "o")
  
  h<-hist(log(GEAs$LASSO_LASSO1se$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue", type = "o")
  
  h<-hist(log(GEAs$LASSO_LM$REACTOME), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green", type = "o")
  
  
  
  
  
  plot(1, type="n", xlab="-log(p-value)", ylab="cummulative counts of number of enriched gene sets", xlim=c(-30, -1), ylim=c(0, 1500), main="GENESIGDB")
  
  h<-hist(log(GEAs$VBSR$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black")
  
  h<-hist(log(GEAs$LASSOmin$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red")
  
  h<-hist(log(GEAs$LASSO1se$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue")
  
  h<-hist(log(GEAs$LM$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green")
  
  h<-hist(log(GEAs$VBSR_VBSR$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black", type = "b")
  
  h<-hist(log(GEAs$VBSR_LASSOmin$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red", type = "b")
  
  h<-hist(log(GEAs$VBSR_LASSO1se$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue", type = "b")
  
  h<-hist(log(GEAs$VBSR_LM$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green", type = "b")
  
  h<-hist(log(GEAs$LASSO_VBSR$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="black", type = "o")
  
  h<-hist(log(GEAs$LASSO_LASSOmin$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="red", type = "o")
  
  h<-hist(log(GEAs$LASSO_LASSO1se$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="blue", type = "o")
  
  h<-hist(log(GEAs$LASSO_LM$GENESIGDB), plot = FALSE, breaks = 50)
  cdf<-cumsum(h$counts)
  lines(h$mids,cdf, col="green", type = "o")
  
  

  
  
  

  #strsplit(names(GEA_REA_LM), "\\.[^\\.]*$")
  
  
  
}
