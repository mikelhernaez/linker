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
GENESETDB_Collections_GeneSymbol_v11<-readMat("~/linker/data/GENESETDB_Collections_GeneSymbol_v11.mat")

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
camodi_input_data<-readMat("~/linker/data/Tumor_OV50_to_R.mat")
rownames(camodi_input_data$geneExpression.matrix)<-unlist(camodi_input_data$geneNames)

lognorm_est_counts<-camodi_input_data$geneExpression.matrix
target_filtered_idx<-camodi_input_data$geneIdx
regulator_filtered_idx<-camodi_input_data$regulatorIdx
##########################################################

Gene_set_Collections<-pathway_genes[c(3,4,5,12)]
#link_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),

testLinker<-LINKER_run(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx, Gene_set_Collections,
                           link_mode=c("LM"),
                           graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                           module_rep="MEAN",
                           NrModules=100, 
                           corrClustNrIter=50,
                           Nr_bootstraps=10,
                           FDR=0.05,
                           NrCores=30)

testNet<-NET_run(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx, Gene_set_Collections,
                  graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                  FDR=0.05,
                  NrCores=30)

LINKER_plot_res_real_data(testLinker$raw_results, file="plot.pdf")
LINKER_plot_GEAs(testLinker$GEAs)
LINKER_plot_graphs_topology(testLinker$graphs)

LINKER_plot_GEAs_norm_regs(testLinker$GEAs, testLinker$graphs,type="EDGE", max_y = c(0.3,0.6,0.6,2), min_x = c(-20,-50,-40,-50))
LINKER_plot_GEAs_norm_regs(testLinker$GEAs, testLinker$graphs,type="REG", max_y = c(.01,.03,.05,.1), min_x = c(-20,-50,-40,-50))

LINKER_plot_GEAs_norm_sets(testLinker$GEAs, testLinker$graphs,type="EDGE", max_y = c(0.3,0.6,0.6,2), min_x = c(-20,-50,-40,-50))
LINKER_plot_GEAs_norm_sets(testLinker$GEAs, testLinker$graphs,type="REG", max_y = c(.01,.03,.05,.1), min_x = c(-20,-50,-40,-50))

LINKER_plot_GEAs_norm_regsets(testLinker$GEAs, testLinker$graphs,type="EDGE", max_y = c(0.3,0.6,0.6,2), min_x = c(-20,-50,-40,-50))
LINKER_plot_GEAs_norm_regsets(testLinker$GEAs, testLinker$graphs,type="REG", max_y = c(.1,.3,.5,1), min_x = c(-20,-50,-40,-50))
#### NOTE: To add NET results to the plots, the $graph and $GEAs of both LinkerOutput and NetOutput must be concatenated. ####
