library("doParallel")
library("igraph")
library("R.matlab")
library("R.utils")
library("vbsr")
library("Matrix")
library("glmnet")
library("colorspace")

source('~/linker/LINKER.R')

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


Foo_lin<-LINKER_run(lognorm_est_counts, target_filtered_idx[1:100], regulator_filtered_idx[1:10], Gene_set_Collections,
                           link_mode=c("VBSR"),
                           graph_mode=c("VBSR"),
                           module_rep="MEAN",
                           NrModules=20, 
                           corrClustNrIter=5,
                           Nr_bootstraps=2,
                           FDR=0.05,
                           NrCores=20)

Foo<-NET_run(lognorm_est_counts, target_filtered_idx[1:100], regulator_filtered_idx[1:10], Gene_set_Collections,
                  graph_mode=c("VBSR"),
                  FDR=0.05,
                  NrCores=30)

LINKER_plot_res_real_data(Foo_lin$raw_results, file="foo.pdf")
LINKER_plot_GEAs(GEAs)
LINKER_plot_graphs_topology(graphs)
