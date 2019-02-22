LINKER_plot_res_real_data<-function(results, file = "plots_real_data.pdf"){
  
  pdf(file)
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0)) 
  
  RsquareAjusted<-list()
  RsquareAjusted_test<-list()
  RsquareAjusted_consensus<-list()
  VarEx<-list()
  VarEx_test<-list()
  VarEx_consensus<-list()
  
  
  for(idx in names(results)){
    res<-results[[idx]]
    Nr_bootstraps<-length(res$bootstrapResults)
    RsquareAjusted[[idx]]<-list(length = Nr_bootstraps)
    RsquareAjusted_test[[idx]]<-list(length = Nr_bootstraps)
    VarEx[[idx]]<-list(length = Nr_bootstraps)
    VarEx_test[[idx]]<-list(length = Nr_bootstraps)
    
    for(i in 1: Nr_bootstraps){
      
      RsquareAjusted[[idx]][[i]]<-rep(res$bootstrapResults[[i]]$trainingStats[,"RsquareAdjusted"],res$bootstrapResults[[i]]$trainingStats[,"nrGen"])
      VarEx[[idx]][[i]]<-rep(res$bootstrapResults[[i]]$trainingStats[,"condition"],res$bootstrapResults[[i]]$trainingStats[,"nrGen"])
      RsquareAjusted_test[[idx]][[i]]<-rep(res$bootstrapTestStats[[i]][,"RsquareAdjusted"],res$bootstrapResults[[i]]$trainingStats[,"nrGen"])
      VarEx_test[[idx]][[i]]<-rep(res$bootstrapTestStats[[i]][,"condition"],res$bootstrapResults[[i]]$trainingStats[,"nrGen"])
    }
  }
  
  
  h<-hist(unlist(RsquareAjusted[[1]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  plot(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "solid", xlab="Adjusted Rsquared",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  h<-hist(unlist(RsquareAjusted[[2]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
  
  legend("topleft", 1.9, 
         c("LASSO","VBSR","Training", "Test"  ), 
         col = c("blue", "red", "black", "black"),
         text.col = "black", lty = c(1,1, 1, 2),
         bg = "gray90", cex = 2)
  
  for(i in 1: Nr_bootstraps){
    h<-hist(RsquareAjusted[[1]][[i]], plot = FALSE, breaks = "FD")
    cdf<-cumsum(h$counts)/sum(h$counts)
    lines(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "solid", xlab="Adjusted Rsquared")
    h<-hist(RsquareAjusted[[2]][[i]], plot = FALSE, breaks = "FD")
    cdf<-cumsum(h$counts)/sum(h$counts)
    lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
    
    h<-hist(RsquareAjusted_test[[1]][[i]], plot = FALSE, breaks = "FD")
    cdf<-cumsum(h$counts)/sum(h$counts)
    lines(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "dotted", xlab="Adjusted Rsquared")
    h<-hist(RsquareAjusted_test[[2]][[i]], plot = FALSE, breaks = "FD")
    cdf<-cumsum(h$counts)/sum(h$counts)
    lines(h$mids, cdf, type = 'l', col = "red4", lty = "dotted")
  }
  
  
  ##### Variance Explained ########
  h<-hist(unlist(VarEx[[1]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  plot(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "solid", xlab="Percentage of variance explained",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  h<-hist(unlist(VarEx[[2]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
  
  legend("bottomright", 1.9, 
         c("LASSO","VBSR","Training", "Test" ), 
         col = c("blue", "red", "black", "black"),
         text.col = "black", lty = c(1,1, 1, 2),
         bg = "gray90", cex = 2)
  
  for(i in 1: Nr_bootstraps){
    h<-hist(VarEx[[1]][[i]], plot = FALSE, breaks = "FD")
    cdf<-cumsum(h$counts)/sum(h$counts)
    lines(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "dotted", xlab="Adjusted Rsquared")
    h<-hist(VarEx[[2]][[i]], plot = FALSE, breaks = "FD")
    cdf<-cumsum(h$counts)/sum(h$counts)
    lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
    
    h<-hist(VarEx_test[[1]][[i]], plot = FALSE, breaks = "FD")
    cdf<-cumsum(h$counts)/sum(h$counts)
    lines(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "dotted", xlab="Adjusted Rsquared")
    h<-hist(VarEx_test[[2]][[i]], plot = FALSE, breaks = "FD")
    cdf<-cumsum(h$counts)/sum(h$counts)
    lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
    
  }
  dev.off()
}


LINKER_plot_graphs_topology<-function(graph_list)
{
  
  col=c("black","red","blue","green","black","red","blue","green","black","red","blue","green")
  pch=c(0,2,5,6,18,19,20,1,3,4,7,8)
  
  module_modes<-names(graph_list)
  par(mfrow=c(5,4))
  for(mode_idx in 1:length(module_modes)){
    graphs<-graph_list[[mode_idx]]
    for(idx in 1:length(graphs)){
      if(class(graphs[[idx]])=="igraph"){
        t<-table(degree(graphs[[idx]]))
        num_edges<-ecount(graphs[[idx]])
      }
      else{
        t<-table(unlist(lapply(graphs[[idx]], function(x) degree(x))))
        num_edges<-sum(unlist(lapply(graphs[[idx]], function(x) ecount(x))))
      }
      plot(as.numeric(names(t)),as.numeric(t)/sum(as.numeric(t)), col=col[idx], 
           pch = 16,log="xy",xlab='Degree',ylab='P(Degree)', main = c(module_modes[mode_idx],names(graphs)[idx],num_edges),xlim=c(1,1000), ylim=c(1/500000,1))
      grid(col = "lightgray", lty = "dotted",
           lwd = par("lwd"), equilogs = TRUE)
    }
  }
  
  par(mfrow=c(5,4))
  for(mode_idx in 1:length(module_modes)){
    graphs<-graph_list[[mode_idx]]
    for(idx in 1:length(graphs)){
      if(class(graphs[[idx]])=="igraph"){
        t<-table(degree(graphs[[idx]], V(graphs[[idx]])[V(graphs[[idx]])$type==0]))
        num_edges<-ecount(graphs[[idx]])
        
      }
      else{
        t<-table(unlist(unname(lapply(graphs[[idx]], function(x) degree(x, V(x)[V(x)$type==0] )))))
        num_edges<-sum(unlist(lapply(graphs[[idx]], function(x) ecount(x))))
        
      }
      plot(as.numeric(names(t)),as.numeric(t)/sum(as.numeric(t)), col=col[idx], 
           pch = 16,log="xy",xlab='Degree',ylab='P(Degree)', main = c(module_modes[mode_idx],names(graphs)[idx],num_edges),xlim=c(1,500), ylim=c(1/500000,1))
      grid(col = "lightgray", lty = "dotted",
           lwd = par("lwd"), equilogs = TRUE)
    }
  }
  
  par(mfrow=c(5,4))
  for(mode_idx in 1:length(module_modes)){
    graphs<-graph_list[[mode_idx]]
    for(idx in 1:length(graphs)){
      if(class(graphs[[idx]])=="igraph"){
        t<-table(degree(graphs[[idx]], V(graphs[[idx]])[V(graphs[[idx]])$type==1]))
        num_edges<-ecount(graphs[[idx]])
        
      }
      else{
        t<-table(unlist(unname(lapply(graphs[[idx]], function(x) degree(x, V(x)[V(x)$type==1] )))))
        num_edges<-sum(unlist(lapply(graphs[[idx]], function(x) ecount(x))))
        
      }
      plot(as.numeric(names(t)),as.numeric(t)/sum(as.numeric(t)), col=col[idx], 
           pch = 16,log="xy",xlab='Degree',ylab='P(Degree)', main = c(module_modes[mode_idx],names(graphs)[idx],num_edges),xlim=c(1,1000), ylim=c(1/500000,1))
      grid(col = "lightgray", lty = "dotted",
           lwd = par("lwd"), equilogs = TRUE)
    }
  }
  
}


LINKER_plot_GEAs<-function(GEAs)
{  
  
  plot(1, type="n", xlab="-log(p-value)", ylab="cummulative counts of number of enriched gene sets", xlim=c(-20, -1), ylim=c(0, 1500), main="BIOCARTA")
  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  
  for(i in 1:length(link_mode)){
    net_mode<-names(GEAs[[ link_mode[i] ]])
    for(j in 1:length(net_mode)){
      h<-hist(log(unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]]$BIOCARTA)), plot = FALSE, breaks = 50)
      cdf<-cumsum(h$counts)
      lines(h$mids,cdf, col=col[i])
      points(h$mids,cdf, col=col[i], pch=j)
    }
  }
  legend("topleft", 1.9, 
         c(link_mode,net_mode), 
         col = c(col,rep.int("black",length(net_mode))),
         text.col = "black", 
         lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
         pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
         bg = "gray90",cex = 2)
  
  
  
  plot(1, type="n", xlab="-log(p-value)", ylab="cummulative counts of number of enriched gene sets", xlim=c(-40, -1), ylim=c(0, 1000), main="KEGG")
  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  
  for(i in 1:length(link_mode)){
    net_mode<-names(GEAs[[ link_mode[i] ]])
    for(j in 1:length(net_mode)){
      h<-hist(log(GEAs[[ link_mode[i] ]][[ net_mode[j] ]]$KEGG), plot = FALSE, breaks = 50)
      cdf<-cumsum(h$counts)
      lines(h$mids,cdf, col=col[i])
      points(h$mids,cdf, col=col[i], pch=j)
    }
  }
  legend("topleft", 1.9, 
         c(link_mode,net_mode), 
         col = c(col,rep.int("black",length(net_mode))),
         text.col = "black", 
         lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
         pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
         bg = "gray90",cex = 2)
  
  
  plot(1, type="n", xlab="-log(p-value)", ylab="cummulative counts of number of enriched gene sets", xlim=c(-40, -1), ylim=c(0, 3000), main="REACTOME")
  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  
  for(i in 1:length(link_mode)){
    net_mode<-names(GEAs[[ link_mode[i] ]])
    for(j in 1:length(net_mode)){
      h<-hist(log(GEAs[[ link_mode[i] ]][[ net_mode[j] ]]$REACTOME), plot = FALSE, breaks = 50)
      cdf<-cumsum(h$counts)
      lines(h$mids,cdf, col=col[i])
      points(h$mids,cdf, col=col[i], pch=j)
    }
  }
  legend("topleft", 1.9, 
         c(link_mode,net_mode), 
         col = c(col,rep.int("black",length(net_mode))),
         text.col = "black", 
         lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
         pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
         bg = "gray90",cex = 2)
  
  
  
  
  
  plot(1, type="n", xlab="-log(p-value)", ylab="cummulative counts of number of enriched gene sets", xlim=c(-40, -1), ylim=c(0, 15000), main="GENESIGDB")
  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  
  for(i in 1:length(link_mode)){
    net_mode<-names(GEAs[[ link_mode[i] ]])
    for(j in 1:length(net_mode)){
      h<-hist(log(GEAs[[ link_mode[i] ]][[ net_mode[j] ]]$GENESIGDB), plot = FALSE, breaks = 50)
      cdf<-cumsum(h$counts)
      lines(h$mids,cdf, col=col[i])
      points(h$mids,cdf, col=col[i], pch=j)
    }
  }
  legend("topleft", 1.9, 
         c(link_mode,net_mode), 
         col = c(col,rep.int("black",length(net_mode))),
         text.col = "black", 
         lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
         pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
         bg = "gray90",cex = 2)
  
  
  
  #strsplit(names(GEA_REA_LM), "\\.[^\\.]*$")
  
}


LINKER_plot_GEAs_boots<-function(GEAs, modules)
{  
  
  plot(1, xlab="Bootstrap Num.", ylab="cummulative counts of new enriched gene sets", xlim=c(0, 3), ylim=c(0, 1500), main="BIOCARTA")
  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  
  for(i in 1:length(link_mode)){
    boot_IDs<-sapply(modules[[ link_mode[i] ]], function(x) x$bootstrap_idx)
    NrBootstraps<-max(boot_IDs)
    net_mode<-names(GEAs[[ link_mode[i] ]])
    for(j in 1:length(net_mode)){
      new_sets<-list()
      union_sets<-c()
      for(boot_id in 1:NrBootstraps){
        boot_sets<-names(unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]]$BIOCARTA[which(boot_IDs==boot_id)]))
        new_sets[[boot_id]]<-setdiff(boot_sets,union_sets)
        union_sets<-union(union_sets,new_sets[[boot_id]])
      }
      lines(1:NrBootstraps,cumsum(sapply(new_sets, length)), col=col[i])
      points(1:NrBootstraps,cumsum(sapply(new_sets, length)), col=col[i], pch=j)
    }
  }
  legend("topleft", 1.9, 
         c(link_mode,net_mode), 
         col = c(col,rep.int("black",length(net_mode))),
         text.col = "black", 
         lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
         pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
         bg = "gray90",cex = 2)
  
  
  
  plot(1, xlab="Bootstrap Num.", ylab="cummulative counts of new enriched gene sets", xlim=c(0, 3), ylim=c(0, 1500), main="KEGG")
  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  
  for(i in 1:length(link_mode)){
    boot_IDs<-sapply(modules[[ link_mode[i] ]], function(x) x$bootstrap_idx)
    NrBootstraps<-max(boot_IDs)
    net_mode<-names(GEAs[[ link_mode[i] ]])
    for(j in 1:length(net_mode)){
      new_sets<-list()
      union_sets<-c()
      for(boot_id in 1:NrBootstraps){
        boot_sets<-names(unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]]$KEGG[which(boot_IDs==boot_id)]))
        new_sets[[boot_id]]<-setdiff(boot_sets,union_sets)
        union_sets<-union(union_sets,new_sets[[boot_id]])
      }
      lines(1:NrBootstraps,cumsum(sapply(new_sets, length)), col=col[i])
      points(1:NrBootstraps,cumsum(sapply(new_sets, length)), col=col[i], pch=j)
    }
  }
  legend("topleft", 1.9, 
         c(link_mode,net_mode), 
         col = c(col,rep.int("black",length(net_mode))),
         text.col = "black", 
         lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
         pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
         bg = "gray90",cex = 2)
  
  
  plot(1, xlab="Bootstrap Num.", ylab="cummulative counts of new enriched gene sets", xlim=c(0, 3), ylim=c(0, 1500), main="REACTOME")
  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  
  for(i in 1:length(link_mode)){
    boot_IDs<-sapply(modules[[ link_mode[i] ]], function(x) x$bootstrap_idx)
    NrBootstraps<-max(boot_IDs)
    net_mode<-names(GEAs[[ link_mode[i] ]])
    for(j in 1:length(net_mode)){
      new_sets<-list()
      union_sets<-c()
      for(boot_id in 1:NrBootstraps){
        boot_sets<-names(unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]]$REACTOME[which(boot_IDs==boot_id)]))
        new_sets[[boot_id]]<-setdiff(boot_sets,union_sets)
        union_sets<-union(union_sets,new_sets[[boot_id]])
      }
      lines(1:NrBootstraps,cumsum(sapply(new_sets, length)), col=col[i])
      points(1:NrBootstraps,cumsum(sapply(new_sets, length)), col=col[i], pch=j)
    }
  }
  legend("topleft", 1.9, 
         c(link_mode,net_mode), 
         col = c(col,rep.int("black",length(net_mode))),
         text.col = "black", 
         lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
         pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
         bg = "gray90",cex = 2)
  
  
  
  
  plot(1, xlab="Bootstrap Num.", ylab="cummulative counts of new enriched gene sets", xlim=c(0, 3), ylim=c(0, 5000), main="GENESIGDB")
  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  
  for(i in 1:length(link_mode)){
    boot_IDs<-sapply(modules[[ link_mode[i] ]], function(x) x$bootstrap_idx)
    NrBootstraps<-max(boot_IDs)
    net_mode<-names(GEAs[[ link_mode[i] ]])
    for(j in 1:length(net_mode)){
      new_sets<-list()
      union_sets<-c()
      for(boot_id in 1:NrBootstraps){
        boot_sets<-names(unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]]$GENESIGDB[which(boot_IDs==boot_id)]))
        new_sets[[boot_id]]<-setdiff(boot_sets,union_sets)
        union_sets<-union(union_sets,new_sets[[boot_id]])
      }
      lines(1:NrBootstraps,cumsum(sapply(new_sets, length)), col=col[i])
      points(1:NrBootstraps,cumsum(sapply(new_sets, length)), col=col[i], pch=j)
    }
  }
  legend("topleft", 1.9, 
         c(link_mode,net_mode), 
         col = c(col,rep.int("black",length(net_mode))),
         text.col = "black", 
         lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
         pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
         bg = "gray90",cex = 2)
  
  
  #strsplit(names(GEA_REA_LM), "\\.[^\\.]*$")
  
}

LINKER_plot_GEAs_norm_regsets<-function(GEAs, graphs, type="EDGE", max_y=c(1,1,1,1), min_x=c(-20,-40,-40,-40), FDR=0.05)
{  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  AUC<-list()
  DB=c("BIOCARTA", "KEGG", "REACTOME", "GENESIGDB")
  for(k in 1:length(DB)){
    AUC[[ DB[k] ]]<-list()
    if(type=="EDGE"){
      plot(1, type="n", xlab="-log(p-value)", ylab="cumm. counts of GEAregset elements per 1K graph edges", xlim=c(min_x[k], -2), ylim=c(0, max_y[k]), main=DB[k])
    }
    else if(type=="REG"){
      plot(1, type="n", xlab="-log(p-value)", ylab="cumm. counts of GEAregset elements per graph regulator", xlim=c(min_x[k], -2), ylim=c(0, max_y[k]), main=DB[k])
    }
    else if(type=="NODE"){
      plot(1, type="n", xlab="-log(p-value)", ylab="cumm. counts of GEAregset elements per 1K graph nodes", xlim=c(min_x[k], -2), ylim=c(0, max_y[k]), main=DB[k])
    }
    else{
      warning("Normalization type not supported")
    }
    
    
    for(i in 1:length(link_mode)){
      net_mode<-names(GEAs[[ link_mode[i] ]])
      AUC[[ DB[k] ]][[ link_mode[i] ]]<-list()
      for(j in 1:length(net_mode)){
        if(length(unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]][[ DB[k] ]]))==0){
          next;
        }
        
        #Normalization values
        if(link_mode[i]=="NET"){
          n_tests<-length(V(graphs[[ link_mode[i] ]][[ net_mode[j] ]])[which(V(graphs[[ link_mode[i] ]][[ net_mode[j] ]])$type==1)])
          if(type=="EDGE"){
            n<-ecount(graphs[[ link_mode[i] ]][[ net_mode[j] ]])/1000
          }
          else if(type=="REG"){
            n<-n_tests
          }
          else if(type=="NODE"){
            n<-vcount(graphs[[ link_mode[i] ]][[ net_mode[j] ]])/1000
          }
          ng<-10
        }
        else{
          n_tests<-sum(unlist(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)length(V(x)[which(V(x)$type==1)]))))
          if(type=="EDGE"){
            n<-sum(unlist(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)ecount(x))))/1000
          }
          else if(type=="REG"){
            n<-sum(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)length(V(x)[which(V(x)$type==1)])))
          }
          else if(type=="NODE"){
            n<-sum(unlist(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)vcount(x))))/1000
          }
          ng<-1
        }
        #print(paste0(link_mode[i]," ",net_mode[j],'         ',ng*n))
        
        #Enrichment tests for all graphs
        GEA_superset<-unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]][[ DB[k] ]])
        GEA_unique<-by(GEA_superset, INDICES=names(GEA_superset), FUN=min)
        
        ##### GEAregset #######
        h<-hist(log(GEA_unique), plot = FALSE, breaks = 50)
        cdf<-cumsum(h$counts)
        
        lines(h$mids + log( n_tests),cdf/n/ng, col=col[i])
        points(h$mids + log( n_tests),cdf/n/ng, col=col[i], pch=j)
        legend("topleft", 1.9, 
               c(link_mode,net_mode), 
               col = c(col,rep.int("black",length(net_mode))),
               text.col = "black", 
               lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
               pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
               bg = "gray90",cex = 2)
        
        q_values<-h$mids + log(n_tests)
        cdf_norm<-cdf[which(q_values<10^FDR)]/n/ng
        AUC[[ DB[k] ]][[ link_mode[i] ]][[ net_mode[j] ]]<-sum(cdf_norm)
      }
    }
  }
  AUC_df<-as.data.frame(sapply(AUC,unlist))
  link_mode<-unlist(strsplit(rownames(AUC_df), "\\.[^\\.]*$"))
  net_mode<-sapply( strsplit(rownames(AUC_df), "\\."), function(x)x[2])
  AUC_df<-cbind(AUC_df, link_mode, net_mode)
  AUC_df<-gather(AUC_df, "DB","auc", c("BIOCARTA", "KEGG", "REACTOME", "GENESIGDB"))
  AUC_df$Rank <- ave( -AUC_df$auc, AUC_df$DB, FUN=rank )
  return(AUC_df)
}

LINKER_plot_GEAs_norm_regs<-function(GEAs, graphs, type="EDGE", max_y=c(1,1,1,1), min_x=c(-20,-40,-40,-40), FDR=0.05)
{  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  DB=c("BIOCARTA", "KEGG", "REACTOME", "GENESIGDB")
  AUC<-list()
  for(k in 1:length(DB)){
    
    AUC[[ DB[k] ]]<-list()
    if(type=="EDGE"){
      plot(1, type="n", xlab="-log(p-value)", ylab="cumm. counts of GEAreg elements per 1K graph edges", xlim=c(min_x[k], -2), ylim=c(0, max_y[k]), main=DB[k])
    }
    else if(type=="REG"){
      plot(1, type="n", xlab="-log(p-value)", ylab="cumm. counts of GEAreg elements per graph regulator", xlim=c(min_x[k], -2), ylim=c(0, max_y[k]), main=DB[k])
    }
    else if(type=="NODE"){
      plot(1, type="n", xlab="-log(p-value)", ylab="cumm. counts of GEAreg elements per 1K graph nodes", xlim=c(min_x[k], -2), ylim=c(0, max_y[k]), main=DB[k])
    }
    else{
      warning("Normalization type not supported")
    }
    
    
    for(i in 1:length(link_mode)){
      net_mode<-names(GEAs[[ link_mode[i] ]])
      AUC[[ DB[k] ]][[ link_mode[i] ]]<-list()
      for(j in 1:length(net_mode)){
        if(length(unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]][[ DB[k] ]]))==0){
          next;
        }
        
        #Normalization values
        if(link_mode[i]=="NET"){
          n_tests<-length(V(graphs[[ link_mode[i] ]][[ net_mode[j] ]])[which(V(graphs[[ link_mode[i] ]][[ net_mode[j] ]])$type==1)])
          if(type=="EDGE"){
            n<-ecount(graphs[[ link_mode[i] ]][[ net_mode[j] ]])/1000
          }
          else if(type=="REG"){
            n<-n_tests
          }
          else if(type=="NODE"){
            n<-vcount(graphs[[ link_mode[i] ]][[ net_mode[j] ]])/1000
          }
          ng<-10
        }
        else{
          n_tests<-sum(unlist(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)length(V(x)[which(V(x)$type==1)]))))
          if(type=="EDGE"){
            n<-sum(unlist(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)ecount(x))))/1000
          }
          else if(type=="REG"){
            n<-sum(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)length(V(x)[which(V(x)$type==1)])))
          }
          else if(type=="NODE"){
            n<-sum(unlist(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)vcount(x))))/1000
          }
          ng<-1
        }
        
        #Enrichment tests for all graphs
        GEA_superset<-unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]][[ DB[k] ]])
        GEA_unique<-by(GEA_superset, INDICES=names(GEA_superset), FUN=min)
        
        ##### GEAreg #######
        GEA_unique_regs<-as.vector(GEA_unique)
        names(GEA_unique_regs)<-unlist(strsplit(names(GEA_unique), "\\.[^\\.]*$"))
        GEA_regs<-by(GEA_unique_regs, INDICES=names(GEA_unique_regs), FUN=min)
        
        h<-hist(log(GEA_regs), plot = FALSE, breaks = 50)
        cdf<-cumsum(h$counts)
        
        lines(h$mids + log( n_tests),cdf/n/ng, col=col[i])
        points(h$mids + log( n_tests),cdf/n/ng, col=col[i], pch=j)
        legend("topleft", 1.9, 
               c(link_mode,net_mode), 
               col = c(col,rep.int("black",length(net_mode))),
               text.col = "black", 
               lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
               pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
               bg = "gray90",cex = 2)
        q_values<-h$mids + log(n_tests)
        cdf_norm<-cdf[which(q_values<10^FDR)]/n/ng
        AUC[[ DB[k] ]][[ link_mode[i] ]][[ net_mode[j] ]]<-sum(cdf_norm)
      }
    }
  }
  AUC_df<-as.data.frame(sapply(AUC,unlist))
  link_mode<-unlist(strsplit(rownames(AUC_df), "\\.[^\\.]*$"))
  net_mode<-sapply( strsplit(rownames(AUC_df), "\\."), function(x)x[2])
  AUC_df<-cbind(AUC_df, link_mode, net_mode)
  AUC_df<-gather(AUC_df, "DB","auc", c("BIOCARTA", "KEGG", "REACTOME", "GENESIGDB"))
  AUC_df$Rank <- ave( -AUC_df$auc, AUC_df$DB, FUN=rank )
  return(AUC_df)
}

LINKER_plot_GEAs_norm_sets<-function(GEAs, graphs, type="EDGE", max_y=c(1,1,1,1), min_x=c(-20,-40,-40,-40), FDR=0.05)
{  
  link_mode<-names(GEAs)
  col=rainbow_hcl(length(link_mode))
  DB=c("BIOCARTA", "KEGG", "REACTOME", "GENESIGDB")
  AUC<-list()
  for(k in 1:length(DB)){
    AUC[[ DB[k] ]]<-list()
    if(type=="EDGE"){
      plot(1, type="n", xlab="-log(p-value)", ylab="cumm. counts of GEAset elements per 1K graph edges", xlim=c(min_x[k], -2), ylim=c(0, max_y[k]), main=DB[k])
    }
    else if(type=="REG"){
      plot(1, type="n", xlab="-log(p-value)", ylab="cumm. counts of GEAset elements per graph regulator", xlim=c(min_x[k], -2), ylim=c(0, max_y[k]), main=DB[k])
    }
    else if(type=="NODE"){
      plot(1, type="n", xlab="-log(p-value)", ylab="cumm. counts of GEAset elements per 1K graph nodes", xlim=c(min_x[k], -2), ylim=c(0, max_y[k]), main=DB[k])
    }
    else{
      warning("Normalization type not supported")
    }
    
    
    for(i in 1:length(link_mode)){
      AUC[[ DB[k] ]][[ link_mode[i] ]]<-list()
      net_mode<-names(GEAs[[ link_mode[i] ]])
      for(j in 1:length(net_mode)){
        if(length(unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]][[ DB[k] ]]))==0){
          next;
        }
        
        #Normalization values
        if(link_mode[i]=="NET"){
          n_tests<-length(V(graphs[[ link_mode[i] ]][[ net_mode[j] ]])[which(V(graphs[[ link_mode[i] ]][[ net_mode[j] ]])$type==1)])
          if(type=="EDGE"){
            n<-ecount(graphs[[ link_mode[i] ]][[ net_mode[j] ]])/1000
          }
          else if(type=="REG"){
            n<-n_tests
          }
          else if(type=="NODE"){
            n<-vcount(graphs[[ link_mode[i] ]][[ net_mode[j] ]])/1000
          }
          ng<-10
        }
        else{
          n_tests<-sum(unlist(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)length(V(x)[which(V(x)$type==1)]))))
          if(type=="EDGE"){
            n<-sum(unlist(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)ecount(x))))/1000
          }
          else if(type=="REG"){
            n<-sum(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)length(V(x)[which(V(x)$type==1)])))
          }
          else if(type=="NODE"){
            n<-sum(unlist(sapply(graphs[[ link_mode[i] ]][[ net_mode[j] ]], function(x)vcount(x))))/1000
          }
          ng<-1
        }
        
        #Enrichment tests for all graphs
        GEA_superset<-unlist(GEAs[[ link_mode[i] ]][[ net_mode[j] ]][[ DB[k] ]])
        GEA_unique<-by(GEA_superset, INDICES=names(GEA_superset), FUN=min)
        
        ##### GEAset #######
        GEA_unique_sets<-as.vector(GEA_unique)
        names(GEA_unique_sets)<-sapply( strsplit(names(GEA_unique), "\\."), function(x)x[2])
        GEA_sets<-by(GEA_unique_sets, INDICES=names(GEA_unique_sets), FUN=min)
        
        h<-hist(log(GEA_sets), plot = FALSE, breaks = 50)
        cdf<-cumsum(h$counts)
        
        lines(h$mids + log( n_tests),cdf/n/ng, col=col[i])
        points(h$mids + log( n_tests),cdf/n/ng, col=col[i], pch=j)
        legend("topleft", 1.9, 
               c(link_mode,net_mode), 
               col = c(col,rep.int("black",length(net_mode))),
               text.col = "black", 
               lty = c(rep.int(1,length(link_mode)),rep.int(-1,length(net_mode))), 
               pch = c(rep.int(-1,length(link_mode)),1:length(net_mode)),
               bg = "gray90",cex = 2)
        q_values<-h$mids + log(n_tests)
        cdf_norm<-cdf[which(q_values<10^FDR)]/n/ng
        AUC[[ DB[k] ]][[ link_mode[i] ]][[ net_mode[j] ]]<-sum(cdf_norm)
      }
    }
  }
  AUC_df<-as.data.frame(sapply(AUC,unlist))
  link_mode<-unlist(strsplit(rownames(AUC_df), "\\.[^\\.]*$"))
  net_mode<-sapply( strsplit(rownames(AUC_df), "\\."), function(x)x[2])
  AUC_df<-cbind(AUC_df, link_mode, net_mode)
  AUC_df<-gather(AUC_df, "DB","auc", c("BIOCARTA", "KEGG", "REACTOME", "GENESIGDB"))
  AUC_df$Rank <- ave( -AUC_df$auc, AUC_df$DB, FUN=rank )
  return(AUC_df)
}

layout.by.attr <- function(graph, wc, cluster.strength=1,layout=layout.auto) {  
  g <- graph.edgelist(get.edgelist(graph)) # create a lightweight copy of graph w/o the attributes.
  E(g)$weight <- 1
  
  attr <- cbind(id=1:vcount(g), val=wc)
  g <- g + vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=cluster.strength)
  
  l <- layout(g, weights=E(g)$weight)[1:vcount(graph),]
  return(l)
}