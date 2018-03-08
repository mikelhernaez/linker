library("doParallel")
library("mvtnorm")
library("igraph")
library("R.matlab")
library("R.utils")
library("vbsr")
library("Matrix")
library("glmnet")

source('~/LINKER_functions.R')
source('~/run_amaretto_lincRNAs.R')

########## Load the gene pathways for the enrichment analysis ##############################
GENESETDB_Collections_GeneSymbol_v11<-readMat("./GENESETDB_Collections_GeneSymbol_v11.mat")

pathway_genes<-list()
for(struct_idx in 1:length(GENESETDB_Collections_GeneSymbol_v11$GENESETDB.Collections))
{
  struct_matrix<-GENESETDB_Collections_GeneSymbol_v11$GENESETDB.Collections[[struct_idx]][[1]][[1]]
  struct_geneList<-GENESETDB_Collections_GeneSymbol_v11$GENESETDB.Collections[[struct_idx]][[1]][[2]]
  struct_pathways<-GENESETDB_Collections_GeneSymbol_v11$GENESETDB.Collections[[struct_idx]][[1]][[3]]
  
  pathway_genes[[struct_idx]]<-list()
  for(geneSet_idx in 1:length(struct_pathways)){
    pathway_genes[[struct_idx]][[struct_pathways[[geneSet_idx]][[1]][1]]]<-sapply(which(struct_matrix[,geneSet_idx]==1), function(x) struct_geneList[[x]][[1]][1])
  }
}
rm(struct_matrix, struct_geneList, struct_pathways, GENESETDB_Collections_GeneSymbol_v11)
###########################################################################################


########## Load the Sleuth data #########
load("/data/olivier/OV.rData")
bootstrap_var<-sleuthSummary$bs_summary$sigma_q_sq
lognorm_est_counts<-sleuthSummary$bs_summary$obs_counts
rm(sleuthSummary)
#########################################


transcript_descr_list<-sapply(names(bootstrap_var), function(x) strsplit(x, "\\|"))
transcript_description<-sapply(transcript_descr_list, function(x) x[[8]])

lincs_idx<-which(transcript_description=="lincRNA")
protein_idx<-which(transcript_description=="protein_coding")

lincs_var<-apply(lognorm_est_counts[lincs_idx,], 1, var)
lincs_bio_var<-lincs_var-bootstrap_var[lincs_idx]

lincs_top_bio_var<-order(lincs_bio_var, decreasing = TRUE)[1:1000]

protein_var<-apply(lognorm_est_counts[protein_idx,], 1, var)
protein_bio_var<-protein_var-bootstrap_var[protein_idx]

protein_top_bio_var<-order(protein_bio_var, decreasing = TRUE)[1:10000]

lincs_filtered_idx <- lincs_idx[lincs_top_bio_var];
protein_filtered_idx <- protein_idx[protein_top_bio_var];

res_OV_100_10<-list()
res_OV_100_10$VBSR<-list()
res_OV_100_10$LASSO<-list()
res_OV_100_10$LM<-list()

for(i in 1:2){
  res_OV_100_10$LASSO[[i]]<-run_linker(lognorm_est_counts, protein_filtered_idx,  lincs_filtered_idx, NrModules=100, module_summary=0,
                          mode="LASSO", used_method="MEAN", corrClustNrIter=100,Nr_bootstraps=10)

  res_OV_100_10$VBSR[[i]]<-run_linker(lognorm_est_counts, protein_filtered_idx,  lincs_filtered_idx, NrModules=100, module_summary=0,
                          mode="VBSR", used_method="MEAN", corrClustNrIter=100,Nr_bootstraps=10)
  
  res_OV_100_10$LM[[i]]<-run_linker(lognorm_est_counts, protein_filtered_idx,  lincs_filtered_idx, NrModules=100, module_summary=0,
                                      mode="LM", used_method="MEAN", corrClustNrIter=100,Nr_bootstraps=10)
}

plot_res_real_data(res_OV_100)

Gene_set_Collections<-pathway_genes[c(3,4,5,12)]
LINKER_plot_enrichement_bootstrap(Gene_set_Collections,res_OV_100[[1]],FDR=0.05)
LINKER_plot_enrichement_bootstrap(Gene_set_Collections,results,FDR=0.05)

enriched_modules_lasso_100<-filter_enriched_modules(Gene_set_Collections,res_OV_100_10$LASSO,FDR=0.05)
enriched_modules_vbsr_100<-filter_enriched_modules(Gene_set_Collections,res_OV_100_10$VBSR,FDR=0.05)
enriched_modules_lm_100<-filter_enriched_modules(Gene_set_Collections,res_OV_100_10$LM,FDR=0.05)

#enriched_graph_lasso_100<-LINKER_createEnrichedGraph(enriched_modules_lasso_100)
#enriched_graph_vbsr_100<-LINKER_createEnrichedGraph(enriched_modules_vbsr_100)


#comm_GEAs_vbsr<-LINKER_compute_enriched_lincs_communities(enriched_modules_vbsr_100, enriched_graph_vbsr_100, lognorm_est_counts)
#comm_GEAs_lasso<-LINKER_compute_enriched_lincs_communities(enriched_modules_lasso_100, enriched_graph_lasso_100, lognorm_est_counts)

graphs<-list()
graphs$VBSR_VBSR<-LINKER_compute_modules_graph(enriched_modules_vbsr_100, lognorm_est_counts, mode="VBSR")
graphs$VBSR_LASSOmin<-LINKER_compute_modules_graph(enriched_modules_vbsr_100, lognorm_est_counts, mode="LASSOmin")
graphs$VBSR_LASSO1se<-LINKER_compute_modules_graph(enriched_modules_vbsr_100, lognorm_est_counts, mode="LASSO1se")
graphs$VBSR_LM<-LINKER_compute_modules_graph(enriched_modules_vbsr_100, lognorm_est_counts, mode="LM")

graphs$LASSO_VBSR<-LINKER_compute_modules_graph(enriched_modules_lasso_100, lognorm_est_counts, mode="VBSR")
graphs$LASSO_LASSOmin<-LINKER_compute_modules_graph(enriched_modules_lasso_100, lognorm_est_counts, mode="LASSOmin")
graphs$LASSO_LASSO1se<-LINKER_compute_modules_graph(enriched_modules_lasso_100, lognorm_est_counts, mode="LASSO1se")
graphs$LASSO_LM<-LINKER_compute_modules_graph(enriched_modules_lasso_100, lognorm_est_counts, mode="LM")

graphs$LM_VBSR<-LINKER_compute_modules_graph(enriched_modules_lm_100, lognorm_est_counts, mode="VBSR")
graphs$LM_LASSOmin<-LINKER_compute_modules_graph(enriched_modules_lm_100, lognorm_est_counts, mode="LASSOmin")
graphs$LM_LASSO1se<-LINKER_compute_modules_graph(enriched_modules_lm_100, lognorm_est_counts, mode="LASSO1se")
graphs$LM_LM<-LINKER_compute_modules_graph(enriched_modules_lm_100, lognorm_est_counts, mode="LM")

graphs$VBSR<-LINKER_compute_graph_all_VBSR(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
graphs$LASSOmin<-LINKER_compute_graph_all_LASSO_min(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
graphs$LASSO1se<-LINKER_compute_graph_all_LASSO_1se(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)
graphs$LM<-LINKER_compute_graph_all_LM(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)

GEAs<-LINKER_compute_graph_enrichment_geneSets_graph_list(pathway_genes,graphs)

LINKER_plot_GEAs(GEAs)
LINKER_plot_graphs_topology(graphs)
  