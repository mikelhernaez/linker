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


########## Load the CaMoDi data #########
camodi_input_data<-readMat("./Tumor_OV50_to_R.mat")
rownames(camodi_input_data$geneExpression.matrix)<-unlist(camodi_input_data$geneNames)

lognorm_est_counts<-camodi_input_data$geneExpression.matrix
protein_filtered_idx<-camodi_input_data$geneIdx
lincs_filtered_idx<-camodi_input_data$regulatorIdx
##########################################################

Gene_set_Collections<-pathway_genes[c(3,4,5,12)]

cl <- makeCluster(16)
registerDoParallel(cl)

camodi_OV<-LINKER_run(lognorm_est_counts, protein_filtered_idx,  lincs_filtered_idx, Gene_set_Collections, NrCores = 30)

NET_networks_creation(lognorm_est_counts, lincs_filtered_idx, protein_filtered_idx)

GEAs<-LINKER_compute_graph_enrichment_geneSets_graph_list(pathway_genes,g)

plot_res_real_data(resC_OV_100_10)
LINKER_plot_GEAs(GEAs)
LINKER_plot_graphs_topology(graphs)
  