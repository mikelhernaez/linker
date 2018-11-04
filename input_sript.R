library("doParallel")
library("igraph")
library("R.matlab")
library("R.utils")
library("vbsr")
library("Matrix")
library("glmnet")
library("colorspace")

source('~/linker/LINKER.R')
source('~/linker/NET_functions.R')
source('~/linker/plot_functions.R')

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


########################################################################

########## Load data #########
camodi_input_data<-readMat("./Tumor_OV50_to_R.mat")
rownames(camodi_input_data$geneExpression.matrix)<-unlist(camodi_input_data$geneNames)

lognorm_est_counts<-camodi_input_data$geneExpression.matrix
target_filtered_idx<-camodi_input_data$geneIdx
regulator_filtered_idx<-camodi_input_data$regulatorIdx
##########################################################

Gene_set_Collections<-pathway_genes[c(3,4,5,12)]


testLinker<-LINKER_run(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx, Gene_set_Collections,
                           link_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                           graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                           module_rep="MEAN",
                           NrModules=100, 
                           corrClustNrIter=50,
                           Nr_bootstraps=10,
                           FDR=0.05,
                           NrCores=20)

testNet<-NET_run(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx, Gene_set_Collections,
                  graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                  FDR=0.05,
                  NrCores=30)

LINKER_plot_res_real_data(testLinker$raw_results, file="plot.pdf")
LINKER_plot_GEAs(testLinker$GEAs)
LINKER_plot_graphs_topology(testLinker$graphs)
#### NOTE: To add NET results to the plots, the $graph and $GEAs of both LinkerOutput and NetOutput must be concatenated. ####
