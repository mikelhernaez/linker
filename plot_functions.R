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
      h<-hist(log(GEAs[[ link_mode[i] ]][[ net_mode[j] ]]$BIOCARTA), plot = FALSE, breaks = 50)
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

layout.by.attr <- function(graph, wc, cluster.strength=1,layout=layout.auto) {  
  g <- graph.edgelist(get.edgelist(graph)) # create a lightweight copy of graph w/o the attributes.
  E(g)$weight <- 1
  
  attr <- cbind(id=1:vcount(g), val=wc)
  g <- g + vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=cluster.strength)
  
  l <- layout(g, weights=E(g)$weight)[1:vcount(graph),]
  return(l)
}