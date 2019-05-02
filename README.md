This is the accompanying repository of the paper:

# Comparison of single and module-based methods for modeling gene regulatory networks
## by Mikel Hernaez, Charles Blatti and Olivier Gevaert


# Abstract

Gene regulatory networks describe the regulatory relationships among genes, and developing methods for reverse engineering these networks is an ongoing challenge in computational biology. The majority of the initially proposed methods for gene regulatory network discovery create a network of genes and then mine it in order to uncover previously unknown regulatory processes. More recent approaches have focused on inferring modules of co-regulated genes, linking these modules with regulator genes and then mining them to discover new molecular biology.

In this work we analyze module-based network approaches to build gene regulatory networks, and compare their performance to the well-established single gene network approaches. In the process, we propose a novel approach to estimate gene regulatory networks drawing from the module-based methods. We show that generating modules of co-expressed genes which are predicted by a sparse set of regulators using a variational bayes method, and then building a bipartite graph on the generated modules using LASSO, yields more informative networks---as measured by the rate of enriched elements and a thorough network topology assessment---than previous single and module-based network approaches. Specifically, the proposed method produces networks closer to a scale-free topology, and the modules show up to 10x more enriched elements than when using single gene networks using TCGA data.

# Methods

The following libraries are needed to run the code:

    library("doParallel")
    library("igraph")
    library("R.matlab")
    library("R.utils")
    library("vbsr")
    library("Matrix")
    library("glmnet")
    library("colorspace")

### Module-based approaches

```LINKER.R``` contain all the functions needed to run the module-based approaches shown in the paper. 

    
The main entry function on ```LINKER.R``` is 

    LINKER_run<-function(
                      lognorm_est_counts, 
                      target_filtered_idx, 
                      regulator_filtered_idx, 
                      Gene_set_Collections,
                      link_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                      graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                      module_rep="MEAN",
                      NrModules=100, 
                      corrClustNrIter=100,
                      Nr_bootstraps=10,
                      FDR=0.05,
                      NrCores=30)
  
Where 

-  ```lognorm_est_counts```: Matrix of log-normalized estimated counts of the gene expression data (Nr Genes x Nr samples)
  
-  ```target_filtered_idx```: Index of the target genes on the ```lognorm_est_counts``` matrix.
-  ```regulator_filtered_idx```: Index of the regulatory genes on the ```lognorm_est_counts``` matrix.
- ```Gene_set_Collections```: Known collection of gene sets for enrichment tests. 
- ```link_mode```: Chosen method(s) to link module eigengenes to regulators. The available options are ```"VBSR", "LASSOmin", "LASSO1se" and "LM"```.
- ```graph_mode```: Chosen method(s) to generate the edges in the bipartite graph. The available options are ```"VBSR", "LASSOmin", "LASSO1se" and "LM"```.
- ```NrModules```: Number of modules that are a priori to be found (note that the final number of modules discovered may differ from this value)  
- ```corrClustNrIter```: Number of iteration for the phase I part of the method.
- ```Nr_bootstraps```: Number of bootstrap of Phase I.
- ```FDR```: The False Discovery Rate correction used for the enrichment analysis.
- ``` NrCores```: Nr of computer cores for the parallel parts of the method. Note that the parallelization is NOT initialized in any of the functions.


### Single gene network approaches

```NET.R``` contain all the functions needed to run the module-based approaches shown in the paper. 

    
The main entry function on ```NET.R``` is 

    NET_run<-function(
                      lognorm_est_counts, 
                      target_filtered_idx, 
                      regulator_filtered_idx, 
                      Gene_set_Collections,
                      graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                      FDR=0.05,
                      NrCores=30)
  
Where 

-  ```lognorm_est_counts```: Matrix of log-normalized estimated counts of the gene expression data (Nr Genes x Nr samples)
  
-  ```target_filtered_idx```: Index of the target genes on the ```lognorm_est_counts``` matrix.
-  ```regulator_filtered_idx```: Index of the regulatory genes on the ```lognorm_est_counts``` matrix.
- ```Gene_set_Collections```: Known collection of gene sets for enrichment tests. 
- ```graph_mode```: Chosen method(s) to generate the edges in the network. The available options are ```"VBSR", "LASSOmin", "LASSO1se" and "LM"```.
- ```FDR```: The False Discovery Rate correction used for the enrichment analysis.
- ``` NrCores```: Nr of computer cores for the parallel parts of the method. Note that the parallelization is NOT initialized in any of the functions.

### Other functions

```plot_functions.R``` contains the needed functions to generate the plots shown in the paper.

```html_functions.R``` contains the needed functions to generate the html summaries of the graph edges.

```simulated_data.R``` contains the functions to generate the simulated data used in the paper and to generate the clustering evaluation plots.


# Data

The data used in the paper can be downloaded from: 

 - Estimated counts and TFs indexes for TCGA OV: [Tumor_HNSC50_to_R.mat](https://github.com/mikelhernaez/linker/blob/master/data/Tumor_HNSC50_to_R.mat)
 - Estimated counts and TFs indexes for TCGA HNSC: [Tumor_OV50_to_R.mat](https://github.com/mikelhernaez/linker/blob/master/data/Tumor_OV50_to_R.mat)
 - Known Gene Sets: [GENESETDB_Collections_GeneSymbol_v11.mat](https://github.com/mikelhernaez/linker/blob/master/data/GENESETDB_Collections_GeneSymbol_v11.mat)
 - TCGA html summaries: [OV](http://donostia.csl.illinois.edu/linker_TCGA/html_Tumor_OV50.tar8855_reg638/index.Tumor_OV50.tar8855_reg638.html) and [HNSC](http://donostia.csl.illinois.edu/linker_TCGA/html_Tumor_HNSC50.tar8791_reg702/index.Tumor_HNSC50.tar8791_reg702.html)

# Example

```input_script.R``` contains an example script.

```html_summary.R``` runs the `create_html_summary()` script on the sample data provided in this repo.
