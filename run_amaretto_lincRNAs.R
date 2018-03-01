run_linker<-function(lognorm_est_counts, protein_filtered_idx, lincs_filtered_idx, NrModules)
{
  
  sample_size<-dim(lognorm_est_counts)[2]
  Nr_bootstraps<-10
  test_size<-round(0.1*sample_size)
  train_size<-round(0.7*sample_size)
  EvaluateTestSet<-list()
  bootstrap_modules<-list()
  bootstrap_results<-list()
  test_samples<-sample(1:sample_size, test_size, replace=F)
  train_val_samples<-setdiff(1:sample_size, test_samples)
  
  Regulator_data_test = t(scale(t(lognorm_est_counts[lincs_filtered_idx,test_samples])))
  MA_matrix_Var_test = t(scale(t(lognorm_est_counts[protein_filtered_idx,test_samples])))
  
  Regulator_data_train_val = t(scale(t(lognorm_est_counts[lincs_filtered_idx,train_val_samples])))
  MA_matrix_Var_train_val = t(scale(t(lognorm_est_counts[protein_filtered_idx,train_val_samples])))
  
  for(boost_idx in 1:Nr_bootstraps)
  {
    
    train_samples<-sample(train_val_samples, train_size, replace=F)
    validation_samples<-setdiff(train_val_samples, train_samples)
    
    Regulator_data_train = t(scale(t(lognorm_est_counts[lincs_filtered_idx,train_samples])))    
    Regulator_data_validation = t(scale(t(lognorm_est_counts[lincs_filtered_idx,validation_samples])))
    
    MA_matrix_Var_train = t(scale(t(lognorm_est_counts[protein_filtered_idx,train_samples])))
    MA_matrix_Var_validation = t(scale(t(lognorm_est_counts[protein_filtered_idx,validation_samples])))
    
  
    # center the data and standardize the regulators
    #MA_matrix_Var = t(scale(t(lognorm_est_counts[protein_filtered_idx,train_samples])))
    #MA_matrix_Var_test = t(scale(t(lognorm_est_counts[protein_filtered_idx,test_samples])))
    
    #MA_matrix_Var = t(apply(lognorm_est_counts[protein_filtered_idx,], 1, function(y) y - mean(y)))
    
    #MA_matrix_Var = scale(lognorm_est_counts[protein_filtered_idx,train_samples])
    #MA_matrix_Var_test = scale(lognorm_est_counts[protein_filtered_idx,test_samples])
  
    LINKERinit<-LINKER_init(MA_matrix_Var = MA_matrix_Var_train, RegulatorData = Regulator_data_train, NrModules = NrModules)
  
    bootstrap_results[[boost_idx]]<-LINKER_corrClust(LINKERinit)
    
    bootstrap_modules[[boost_idx]]<-list()
    for(i in 1:bootstrap_results[[boost_idx]]$NrModules){
      
      module_genes<-bootstrap_results[[boost_idx]]$AllGenes[which(bootstrap_results[[boost_idx]]$ModuleMembership==i)]
      module_gene_description_list<-sapply(module_genes, function(x) strsplit(x, "\\|"))
      bootstrap_modules[[boost_idx]][[i]]<-sapply(module_gene_description_list, function(x) x[[6]])
    }
  
    EvaluateTestSet[[boost_idx]] <- LINKER_EvaluateTestSet(bootstrap_results[[boost_idx]],MA_matrix_Var_validation,Regulator_data_validation)
  
    printf("Bootstrap %d, NrModules %d:\n", boost_idx, bootstrap_results[[boost_idx]]$NrModules)
    
    print(apply(EvaluateTestSet[[boost_idx]], 2, mean))
  }
  
  #g<-create_bootstrap_graph(bootstrap_results, bootstrap_modules)
    
  corrClust_results<-list(bootstrap_modules, EvaluateTestSet, bootstrap_results)
  
  bootstrap_matrix<-LINKER_filterModuleBootstraps(corrClust_results)
  
  consensus_modules<-LINKER_consensusClustering(bootstrap_matrix)
  
  regulators<-LINKER_compute_lincRegulators(MA_matrix_Var_train_val, Regulator_data_train_val, consensus_modules)
  
  testResults<-LINKER_EvaluateTestSet(regulators,MA_matrix_Var_test,Regulator_data_test)
  
  printf("NrModules %d:\n",nrow(testResults))
  print(apply(testResults, 2, mean))
  
  return(list(testResults, regulators))
  
}









