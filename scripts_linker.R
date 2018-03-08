# Get data from CaMoDi

camodi_input_data<-readMat("./Tumor_OV50_to_R.mat")
rownames(camodi_input_data$geneExpression.matrix)<-unlist(camodi_input_data$geneNames)

res<-list()
res[[1]]<-list()
res[[2]]<-list()
sim<-list()
for(i in 1:2){
  
  sim[[i]]<-generate_sim_data(camodi_input_data$geneExpression.matrix[camodi_input_data$regulatorIdx,])
  
  res[[1]][[i]]<-run_linker(sim[[i]][[1]], sim[[i]][[2]], sim[[i]][[3]], NrModules=40, module_summary=1,
                            mode="LASSO", used_method="MEAN", corrClustNrIter=30,Nr_bootstraps=20)
  
  res[[2]][[i]]<-run_linker(sim[[i]][[1]], sim[[i]][[2]], sim[[i]][[3]], NrModules=40, module_summary=1,
                            mode="VBSR", used_method="MEAN", corrClustNrIter=30,Nr_bootstraps=20)
  
  #res[[1]][[i]]<-run_linker(camodi_input_data$geneExpression.matrix, camodi_input_data$geneIdx, camodi_input_data$regulatorIdx, NrModules=100, module_summary=1,
  #                          mode="LASSO", used_method="MEAN", corrClustNrIter=30,Nr_bootstraps=10)
  
  #res[[2]][[i]]<-run_linker(camodi_input_data$geneExpression.matrix, camodi_input_data$geneIdx, camodi_input_data$regulatorIdx, NrModules=100, module_summary=1,
  #                          mode="VBSR", used_method="MEAN", corrClustNrIter=30,Nr_bootstraps=10)
  
}

plot_clust_eval(sim, res, 2)

evaluate_fit(res, sim, used_method = "MEAN")

plot_res_real_data(res)


yplot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", xlab = "x", ylab = "Density")
#modules<-c(25,50,100,150,200)

modules<-c(25, 50, 100)
iters<-c(5, 15, 30, 50)
mode<-c("LASSO", "VBSR")
idx<-1
res_camodi_SVD_ABS_OV_train<-list()
for(idx_modules in modules){
  for(idx_iters in iters){
    for(idx_mode in mode){
      res_camodi_SVD_ABS_OV_train[[idx]]<-run_linker(camodi_input_data$geneExpression.matrix, camodi_input_data$geneIdx, camodi_input_data$regulatorIdx, NrModules=idx_modules, module_summary=0,
                             mode=idx_mode, corrClustNrIter=idx_iters,Nr_bootstraps=10)
      idx<-idx+1
    }
  }
}

plot(range(1, d2$x), range(d1$y, d2$y), type = "n", xlab = "x", ylab = "Density")


d2<-(density(resVBSR200[[1]][,5]))
lines(d2, col = "black")

res<-resVBSR
res<-res[[1]]
ordered_res<-res[order(res[,5]),]
first_module<-length(which(cumsum(ordered_res[,2])<0.2*sum(ordered_res[,2])))
res_filt<-ordered_res[first_module:nrow(ordered_res),]
print(apply(res_filt, 2, mean))



counts<-matrix(0,nrow = 24, ncol = 20-1)
for(j in 1:24){
  for(i in 1:10){
    h<-hist(rep(res_camodi_SVD_ABS_OV[[j]][[i]][,5],res_camodi_SVD_ABS_OV[[j]][[i]][,2]), breaks=seq(-4,1, len=20),plot=FALSE)
    counts[j,]<-counts[j,]+h$counts
  }
  counts[j,]<-counts[j,]/sum(counts[j,])
}
cdfs<-apply(counts,1,cumsum)

cdfs_abs<-cdfs[,c(rbind(seq(3,24,by=4),seq(4,24,by=4)))]

counts<-matrix(0,nrow = 24, ncol = 20-1)
for(j in 1:24){
  for(i in 1:10){
    h<-hist(rep(res_camodi_SVD_COR_OV[[j]][[i]][,5],res_camodi_SVD_COR_OV[[j]][[i]][,2]), breaks=seq(-4,1, len=20),plot=FALSE)
    counts[j,]<-counts[j,]+h$counts
  }
  counts[j,]<-counts[j,]/sum(counts[j,])
}
cdfs<-apply(counts,1,cumsum)

cdfs_cor<-cdfs[,c(rbind(seq(3,24,by=4),seq(4,24,by=4)))]


matplot(h$mids,cbind(cdfs_abs, cdfs_cor), type="b")

