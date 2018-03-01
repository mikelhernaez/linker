verify_new_run<-function(amaretto_results, new_amaretto_results, Module_number)
{
  #### Convert results to modules for enrichment analysis #####
  old_modules<-list()
  for(Module_idx in 1:50){
    Module_protein_coding_genes<-names(which(amaretto_results[[3]]$ModuleMembership[,]==Module_idx))
    Module_protein_coding_gene_list<-sapply(Module_protein_coding_genes, function(x) strsplit(x, "\\|"))
    Module_protein_coding_genes<-sapply(Module_protein_coding_gene_list, function(x) x[[6]])
    Module_protein_coding_genes<-unname(Module_protein_coding_genes)
    old_modules[[Module_idx]]<-Module_protein_coding_genes
  }
  
  #### Compare the new results to the old modules
  total_genes<-10000
  Module<-new_amaretto_results[[1]][[Module_number]]
  i<-1
  under_enrichment_pvalues<-numeric()
  over_enrichment_pvalues<-numeric()
  for(pathway_idx in 1:length(old_modules)){
    
    white_balls<-length(old_modules[[pathway_idx]])
    black_balls<-total_genes - white_balls
    
    drawn<-length(Module)
    drawn_whites<-length(which(old_modules[[pathway_idx]] %in% Module))
    
    over_enrichment_pvalues[i]<- phyper(drawn_whites-1, white_balls, black_balls, drawn, lower.tail = FALSE, log.p = FALSE)
    i<-i+1
  }
  
  over_adjp<-p.adjust(over_enrichment_pvalues,'holm')
  
  enriched_modules<-which(over_adjp<1e-10)
  
  return(enriched_modules)
  
}