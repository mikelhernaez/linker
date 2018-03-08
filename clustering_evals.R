evaluate_clusterings<-function(true_assign, bootstrap_results, consensus_results)
  {
  
    cluster_metrics = mat.or.vec(length(bootstrap_results)+1,4)
    boot_assign<-sapply(bootstrap_results, function(x) x$ModuleMembership)
    
    rownames(boot_assign)<-rownames(bootstrap_results[[1]]$ModuleMembership)
  
    for(boot_idx in 1:length(bootstrap_results)){
    
      contingency_matrix<-matrix(0,nrow = length(unique(true_assign)), ncol= max(unique(boot_assign[,boot_idx])))
    
      for(true_C_idx in 1: nrow(contingency_matrix)){
        true_gen<-names(which(true_assign==true_C_idx))
        boot_gen<-which(names(boot_assign[,boot_idx])%in%true_gen)
        count_table<-as.matrix(table(boot_assign[boot_gen,boot_idx]))
        contingency_matrix[true_C_idx,as.numeric(rownames(count_table))]<-count_table
      }
      
      N<-sum(contingency_matrix)
      Pkk<-contingency_matrix/N
      boot_nk<-apply(contingency_matrix, 2, sum)
      boot_Pk<-boot_nk/N
      true_nk<-apply(contingency_matrix, 1, sum)
      true_Pk<-true_nk/N
    
      boot_H<- -sum(sapply(boot_Pk, function(x) {ifelse(x>0,x*log2(x),0)}))
      true_H<- -sum(sapply(true_Pk, function(x) {ifelse(x>0,x*log2(x),0)}))
      
      MI<-0
      for(i in 1:nrow(contingency_matrix)){
        for(j in 1:ncol(contingency_matrix)){
          if(Pkk[i,j]==0){
            next
          }
          MI<-MI+ Pkk[i,j] * log2( Pkk[i,j]/(boot_Pk[j]*true_Pk[i]) )
        }
      }
     
      cluster_metrics[boot_idx,1]<-boot_H + true_H - 2*MI
      cluster_metrics[boot_idx,2]<-MI/sqrt(boot_H*true_H)
      cluster_metrics[boot_idx,3]<-1/(2*N)*(2*N - sum(apply(contingency_matrix, 1, max)) - sum(apply(contingency_matrix, 2, max)))
      
      true_pair<-sum(sapply(true_nk, function(x) choose(x,2)))
      boot_pair<-sum(sapply(boot_nk, function(x) choose(x,2)))
      all_pair<-sum(choose(contingency_matrix, 2))
      
      cluster_metrics[boot_idx,4]<-( all_pair - (true_pair*boot_pair)/choose(N,2) ) / ( 0.5*(true_pair + boot_pair) - ((true_pair*boot_pair)/choose(N,2)) )
      colnames(cluster_metrics)<-c("VI Distance", "Adjusted MI","Van Dogen Distance", "Adjusted Rand Index")
    }
    
    ### Consensus results ###
    consensus_assign<-consensus_results$ModuleMembership
    
    rownames(consensus_results$ModuleMembership)<-rownames(bootstrap_results[[1]]$ModuleMembership)
    
    rownames(consensus_assign)<-rownames(consensus_results$ModuleMembership)
    
    contingency_matrix<-matrix(0,nrow = length(unique(true_assign)), ncol= max(unique(consensus_assign)))
    
    for(true_C_idx in 1: nrow(contingency_matrix)){
      true_gen<-names(which(true_assign==true_C_idx))
      boot_gen<-which(rownames(consensus_assign)%in%true_gen)
      count_table<-as.matrix(table(consensus_assign[boot_gen]))
      contingency_matrix[true_C_idx,as.numeric(rownames(count_table))]<-count_table
    }
    
    N<-sum(contingency_matrix)
    Pkk<-contingency_matrix/N
    boot_nk<-apply(contingency_matrix, 2, sum)
    boot_Pk<-boot_nk/N
    true_nk<-apply(contingency_matrix, 1, sum)
    true_Pk<-true_nk/N
    
    boot_H<- -sum(sapply(boot_Pk, function(x) {ifelse(x>0,x*log2(x),0)}))
    true_H<- -sum(sapply(true_Pk, function(x) {ifelse(x>0,x*log2(x),0)}))
    
    MI<-0
    for(i in 1:nrow(contingency_matrix)){
      for(j in 1:ncol(contingency_matrix)){
        if(Pkk[i,j]==0){
          next
        }
        MI<-MI+ Pkk[i,j] * log2( Pkk[i,j]/(boot_Pk[j]*true_Pk[i]) )
      }
    }
    
    cluster_metrics[boot_idx+1,1]<-boot_H + true_H - 2*MI
    cluster_metrics[boot_idx+1,2]<-MI/sqrt(boot_H*true_H)
    cluster_metrics[boot_idx+1,3]<-1/(2*N)*(2*N - sum(apply(contingency_matrix, 1, max)) - sum(apply(contingency_matrix, 2, max)))
    
    true_pair<-sum(sapply(true_nk, function(x) choose(x,2)))
    boot_pair<-sum(sapply(boot_nk, function(x) choose(x,2)))
    all_pair<-sum(choose(contingency_matrix, 2))
    
    cluster_metrics[boot_idx+1,4]<-( all_pair - (true_pair*boot_pair)/choose(N,2) ) / ( 0.5*(true_pair + boot_pair) - ((true_pair*boot_pair)/choose(N,2)) )
    
    dimnames(cluster_metrics) <- list(rownames(cluster_metrics, do.NULL = FALSE, prefix = "bootstrap_"), c("VI" , "norm_MI", "van_dogen", "ARI"))
    
    
    return(cluster_metrics)
}

plot_clust_eval<-function(sim,res,Num_sims){
  
  
  pdf("clusts.pdf")
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0)) 
  
  
  clust_stats_lasso<-evaluate_clusterings(sim[[1]][[5]], res[[1]][[1]]$bootstrapResults, res[[1]][[1]]$consensusResults)
  clust_stats_vbsr<-evaluate_clusterings(sim[[1]][[5]], res[[2]][[1]]$bootstrapResults, res[[2]][[1]]$consensusResults)
  
  for(i in 2:Num_sims){
    
    clust_stats_lasso<-rbind(clust_stats_lasso,evaluate_clusterings(sim[[i]][[5]], res[[1]][[i]]$bootstrapResults, res[[1]][[i]]$consensusResults))
    clust_stats_vbsr<-rbind(clust_stats_vbsr,evaluate_clusterings(sim[[i]][[5]], res[[2]][[i]]$bootstrapResults, res[[2]][[i]]$consensusResults))
  }
  consensus_stats_lasso<-clust_stats_lasso[which(rownames(clust_stats_lasso)=="bootstrap_21"),]
  consensus_stats_vbsr<-clust_stats_vbsr[which(rownames(clust_stats_vbsr)=="bootstrap_21"),]
  
  clust_stats_lasso<-clust_stats_lasso[-which(rownames(clust_stats_lasso)=="bootstrap_21"),]
  clust_stats_vbsr<-clust_stats_vbsr[-which(rownames(clust_stats_vbsr)=="bootstrap_21"),]
  
  #par(mfrow=c(2,2))
  #pdf("clust_evals.pdf")
  plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(0,2.5),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
       xlab = "", ylab = "Variation of Information",xaxt='n')
  vioplot(clust_stats_lasso[,1], clust_stats_vbsr[,1],col=rgb(0.1,0.4,0.7,0.7),names=c("LASSO","VBSR"),add=TRUE)
  axis(side=1,at=1:2,labels=c("LASSO","VBSR"),cex.lab=2, cex.axis=2)
  #points(rep(1,time=length(consensus_stats_lasso[,1])),consensus_stats_lasso[,1],col = "red")
  #points(rep(2,time=length(consensus_stats_vbsr[,1])),consensus_stats_vbsr[,1],col = "red")
  
  #vioplot(clust_stats_lasso[,2], clust_stats_vbsr[,2],col=rgb(0.1,0.4,0.7,0.7) , names=c("LASSO","VBSR"),
  #        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
  #        ylab =colnames(clust_stats_lasso)[2])
  #points(rep(1,time=length(consensus_stats_lasso[,1])),consensus_stats_lasso[,2],col = "red")
  #points(rep(2,time=length(consensus_stats_vbsr[,1])),consensus_stats_vbsr[,2],col = "red")
  #title(main=colnames(clust_stats_lasso)[2])

  #vioplot(clust_stats_lasso[,3], clust_stats_vbsr[,3],col=rgb(0.1,0.4,0.7,0.7) , names=c("LASSO","VBSR"),
  #        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
  #        ylab =colnames(clust_stats_lasso)[3])
  #points(rep(1,time=length(consensus_stats_lasso[,1])),consensus_stats_lasso[,3],col = "red")
  #points(rep(2,time=length(consensus_stats_vbsr[,1])),consensus_stats_vbsr[,3],col = "red")
  #title(main=colnames(clust_stats_lasso)[3])
  
  plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(0.5,1),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
       xlab = "", ylab = "Adjusted Rand Index",xaxt='n')
  vioplot(clust_stats_lasso[,4], clust_stats_vbsr[,4],col=rgb(0.1,0.4,0.7,0.7),names=c("LASSO","VBSR"),add=TRUE)
  axis(side=1,at=1:2,labels=c("LASSO","VBSR"),cex.lab=2, cex.axis=2)
  
  #vioplot(clust_stats_lasso[,4], clust_stats_vbsr[,4],col=rgb(0.1,0.4,0.7,0.7) , names=c("LASSO","VBSR"),
  #        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
  #        ylab =colnames(clust_stats_lasso)[4])
  #points(rep(1,time=length(consensus_stats_lasso[,1])),consensus_stats_lasso[,4],col = "red")
  #points(rep(2,time=length(consensus_stats_vbsr[,1])),consensus_stats_vbsr[,4],col = "red")
  #title(main=colnames(clust_stats_lasso)[4])
  dev.off()
}

evaluate_fit<-function(res, sim, used_method = "LINKER"){
  
  NrSims<-length(sim)
  
  pdf("plots.pdf")
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0)) 
  
  sim_results<-evaluate_sim_data(sim, used_method = used_method)
  true_conf_matrix<-sim_results$conf_matrix
  
  true_reg<-apply(sim[[1]][[4]], 1, function(x)sum(x!=0))
  
  Nr_bootstraps<-length(res[[1]][[1]]$bootstrapResults)
  confussion_matrices<-list()
  confussion_matrices_consensus<-list()
  sim_betas<-list()
  true_betas<-list()
  sim_betas_consensus<-list()
  true_betas_consensus<-list()
  RsquareAjusted<-list()
  RsquareAjusted_consensus<-list()
  VarEx<-list()
  VarEx_consensus<-list()
  
  
  
  for(idx in 1:2){
    matrix_idx<-1
    confussion_matrices[[idx]]<-matrix(nrow = 4, ncol = Nr_bootstraps*NrSims)
    rownames(confussion_matrices[[idx]])<-c("TP", "FP", "TN", "FN")
    confussion_matrices_consensus[[idx]]<-matrix(nrow = 4, ncol = NrSims)
    rownames(confussion_matrices_consensus[[idx]])<-c("TP", "FP", "TN", "FN")
    sim_betas[[idx]]<-list(length = NrSims*Nr_bootstraps)
    true_betas[[idx]]<-list(length = NrSims*Nr_bootstraps)
    sim_betas_consensus[[idx]]<-list(length = NrSims)
    true_betas_consensus[[idx]]<-list(length = NrSims)
    RsquareAjusted[[idx]]<-list(length = NrSims*Nr_bootstraps)
    RsquareAjusted_consensus[[idx]]<-list(length = NrSims)
    VarEx[[idx]]<-list(length = NrSims*Nr_bootstraps)
    VarEx_consensus[[idx]]<-list(length = NrSims)
    
    for(j in 1:NrSims){
      true_reg<-apply(sim[[j]][[4]], 1, function(x)sum(x!=0))
      
      for(i in 1: Nr_bootstraps){
        sim_reg<-apply(res[[idx]][[j]]$bootstrapResults[[i]]$RegulatoryPrograms, 2, function(x)sum(x!=0))
        TN=sum( (sim_reg == true_reg) & (true_reg==0) )
        TP=sum( ((sim_reg == true_reg) & (true_reg!=0)) * true_reg)
        FP=sum( ((sim_reg != true_reg) & (true_reg==0)) * sim_reg)
        FN=sum( ((sim_reg != true_reg) & (sim_reg==0)) * true_reg)
        
        # We need to adjust for the mistmatched ones
        reminder_mm<-abs(sim_reg-true_reg)
        TP = TP + sum( ((sim_reg != true_reg) & (sim_reg!=0) & (true_reg!=0)) * pmin(sim_reg,true_reg) )
        FP = FP + sum( ((sim_reg > true_reg) & (sim_reg!=0) & (true_reg!=0)) * reminder_mm )
        FN = FN + sum( ((sim_reg < true_reg) & (sim_reg!=0) & (true_reg!=0)) * reminder_mm )
        
        confussion_matrices[[idx]][1,matrix_idx]<-TP
        confussion_matrices[[idx]][2,matrix_idx]<-FP
        confussion_matrices[[idx]][3,matrix_idx]<-TN
        confussion_matrices[[idx]][4,matrix_idx]<-FN

        
        single_tp_idx<-which((sim_reg == true_reg) & (true_reg==1))
        sim_betas[[idx]][[matrix_idx]]<-scale(apply(res[[idx]][[j]]$bootstrapResults[[i]]$RegulatoryPrograms[,single_tp_idx], 2, sum), center=FALSE)
        true_betas[[idx]][[matrix_idx]]<-scale(apply(sim[[j]][[4]][single_tp_idx,], 1, sum), center=FALSE)

        RsquareAjusted[[idx]][[matrix_idx]]<-rep(res[[idx]][[j]]$bootstrapResults[[i]]$trainingStats[,"RsquareAdjusted"],res[[idx]][[j]]$bootstrapResults[[i]]$trainingStats[,"nrGen"])
        VarEx[[idx]][[matrix_idx]]<-rep(res[[idx]][[j]]$bootstrapResults[[i]]$trainingStats[,"condition"],res[[idx]][[j]]$bootstrapResults[[i]]$trainingStats[,"nrGen"])
        
        
        matrix_idx<-matrix_idx + 1
        
        # Data<-sim[[j]][[1]][sim[[j]][[2]] , ]
        # RegulatorData<-sim[[j]][[1]][sim[[j]][[3]] , ]
        # Clusters<-res[[idx]][[j]]$bootstrapResults[[i]]$ModuleMembership
        # ClusterIDs<-unique(Clusters)
        # NrClusters<-length(ClusterIDs)
        # 
        # for (clustIdx in 1:NrClusters){
        # 
        #   CurrentClusterPositions = which(Clusters==ClusterIDs[clustIdx])
        #   nrGenesInClusters = length(CurrentClusterPositions)
        # 
        #   if(used_method=="LINKER"){
        #     gene_svd<-svd(Data[CurrentClusterPositions,])
        #     y<-gene_svd$v[,1]
        #     if(cor(apply(Data[CurrentClusterPositions,],2,mean), y)<0){
        #       y<- -y
        #     }
        #   }
        #   else{
        #     y<-apply(Data[CurrentClusterPositions,],2,mean)
        #   }
        # 
        #   RsquareAjusted[[idx]][clustIdx]<-compute_adj_r2(betas=res[[idx]][[j]]$bootstrapResults[[i]]$RegulatoryPrograms[ClusterIDs[clustIdx],], y, regs=RegulatorData)
        # }
        
      }
      
      sim_reg<-apply(res[[idx]][[j]]$consensusResults$RegulatoryPrograms, 2, function(x)sum(x!=0))
      TN=sum( (sim_reg == true_reg) & (true_reg==0) )
      TP=sum( ((sim_reg == true_reg) & (true_reg!=0)) * true_reg)
      FP=sum( ((sim_reg != true_reg) & (true_reg==0)) * sim_reg)
      FN=sum( ((sim_reg != true_reg) & (sim_reg==0)) * true_reg)
      
      # We need to adjust for the mistmatched ones
      reminder_mm<-abs(sim_reg-true_reg)
      TP = TP + sum( ((sim_reg != true_reg) & (sim_reg!=0) & (true_reg!=0)) * pmin(sim_reg,true_reg) )
      FP = FP + sum( ((sim_reg > true_reg) & (sim_reg!=0) & (true_reg!=0)) * reminder_mm )
      FN = FN + sum( ((sim_reg < true_reg) & (sim_reg!=0) & (true_reg!=0)) * reminder_mm )
      
      confussion_matrices_consensus[[idx]][1,j]<-TP
      confussion_matrices_consensus[[idx]][2,j]<-FP
      confussion_matrices_consensus[[idx]][3,j]<-TN
      confussion_matrices_consensus[[idx]][4,j]<-FN
      
      single_tp_idx<-which((sim_reg == true_reg) & (true_reg==1))
      sim_betas_consensus[[idx]][[j]]<-scale(apply(res[[idx]][[j]]$consensusResults$RegulatoryPrograms[,single_tp_idx], 2, sum), center=FALSE)
      true_betas_consensus[[idx]][[j]]<-scale(apply(sim[[j]][[4]][single_tp_idx,], 1, sum), center=FALSE)
      
      RsquareAjusted_consensus[[idx]][[j]]<-rep(res[[idx]][[j]]$consensusTestStats[,"RsquareAdjusted"],res[[idx]][[j]]$consensusTestStats[,"nrGen"])
      VarEx_consensus[[idx]][[j]]<-rep(res[[idx]][[j]]$consensusTestStats[,"condition"],res[[idx]][[j]]$consensusTestStats[,"nrGen"])
      
    }
    
  }
  
  sim_betas_lasso<-unlist(sim_betas[[1]])
  sim_betas_lasso_cons<-unlist(sim_betas_consensus[[1]])
  sim_betas_vbsr<-unlist(sim_betas[[2]])
  sim_betas_vbsr_cons<-unlist(sim_betas_consensus[[2]])
  true_reg_lasso<-unlist(true_betas[[1]])
  true_reg_vbsr<-unlist(true_betas[[2]])
  true_reg_lasso_cons<-unlist(true_betas_consensus[[1]])
  true_reg_vbsr_cons<-unlist(true_betas_consensus[[2]])
  

  
  model_all<-lm(true_reg_lasso ~ sim_betas_lasso)
  plot(sim_betas_lasso,true_reg_lasso, pch=20, col="blue",
       xlab= expression(paste("true weights ", beta)), ylab=expression(paste("estimated weights ", beta)),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2);
  abline(model_all$coefficients, col="darkblue");
  #model_all<-lm(true_reg_lasso_cons ~ sim_betas_lasso_cons)
  #points(sim_betas_lasso_cons,true_reg_lasso_cons, pch=20, col="red");abline(model_all$coefficients, col="darkred");
  legend("topleft", 1.9, 
         c(paste0("MSE: ", round(sum((true_reg_lasso-sim_betas_lasso)^2)/length(true_reg_lasso), digits = 5) )), 
         col = c("blue"),
         text.col = "black", lty = c(-1), pch = c(20),
         cex = 2,
        bg = "gray90")
  # legend("topleft", 1.9, 
  #        c(paste0("TP MSE: ",sum((true_reg_lasso-sim_betas_lasso)^2)/length(true_reg_lasso))), 
  #        col = c("blue"),
  #        text.col = "black", lty = c(-1), pch = c(20),
  #        bg = "gray90")
  abline(a=0,b=1, col="black")
  
  model_all<-lm(true_reg_vbsr ~ sim_betas_vbsr)
  plot(sim_betas_vbsr,true_reg_vbsr, pch=20, col="blue",
       xlab= expression(paste("true weights ", beta)), ylab=expression(paste("estimated weights ", beta)),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2);
  abline(model_all$coefficients, col="darkblue");
  #model_all<-lm(true_reg_vbsr_cons ~ sim_betas_vbsr_cons)
  #points(sim_betas_vbsr_cons,true_reg_vbsr_cons, pch=20, col="red");abline(model_all$coefficients, col="darkred");
  # legend("topleft", 1.9, 
  #        c(paste0("TP MSE: ",sum((true_reg_vbsr-sim_betas_vbsr)^2)/length(true_reg_vbsr)),paste0("Cons TP MSE: ",sum((true_reg_vbsr_cons-sim_betas_vbsr_cons)^2)/length(true_reg_vbsr_cons))), 
  #        col = c("blue", "red"),
  #        text.col = "black", lty = c(-1,-1), pch = c(20, 20),
  #        bg = "gray90")
  legend("topleft",1.9,  
         c(paste0("MSE: ", round(sum((true_reg_vbsr-sim_betas_vbsr)^2)/length(true_reg_vbsr),digits=5) )), 
         col = c("blue"),
         text.col = "black", lty = c(-1), pch = c( 20),
         bg = "gray90",
         cex = 2)
  abline(a=0,b=1, col="black")
  
  Prec_L<-confussion_matrices[[1]]["TP",]/(confussion_matrices[[1]]["TP",]+confussion_matrices[[1]]["FP",])
  Prec_V<-confussion_matrices[[2]]["TP",]/(confussion_matrices[[2]]["TP",]+confussion_matrices[[2]]["FP",])
  Prec_L_cons<-confussion_matrices_consensus[[1]]["TP",]/(confussion_matrices_consensus[[1]]["TP",]+confussion_matrices_consensus[[1]]["FP",])
  Prec_V_cons<-confussion_matrices_consensus[[2]]["TP",]/(confussion_matrices_consensus[[2]]["TP",]+confussion_matrices_consensus[[2]]["FP",])
  Prec_L_trueClust<-true_conf_matrix[[1]]["TP",]/(true_conf_matrix[[1]]["TP",]+true_conf_matrix[[1]]["FP",])
  Prec_V_trueClust<-true_conf_matrix[[2]]["TP",]/(true_conf_matrix[[2]]["TP",]+true_conf_matrix[[2]]["FP",])
  # vioplot(Prec_L, Prec_V,col=rgb(0.1,0.4,0.7,0.7) , names=c("LASSO","VBSR"))
  # title(main="Precision")
  # 
  Recall_L<-confussion_matrices[[1]]["TP",]/(confussion_matrices[[1]]["TP",]+confussion_matrices[[1]]["FN",])
  Recall_V<-confussion_matrices[[2]]["TP",]/(confussion_matrices[[2]]["TP",]+confussion_matrices[[2]]["FN",])
  Recall_L_cons<-confussion_matrices_consensus[[1]]["TP",]/(confussion_matrices_consensus[[1]]["TP",]+confussion_matrices_consensus[[1]]["FN",])
  Recall_V_cons<-confussion_matrices_consensus[[2]]["TP",]/(confussion_matrices_consensus[[2]]["TP",]+confussion_matrices_consensus[[2]]["FN",])
  Recall_L_trueClust<-true_conf_matrix[[1]]["TP",]/(true_conf_matrix[[1]]["TP",]+true_conf_matrix[[1]]["FN",])
  Recall_V_trueClust<-true_conf_matrix[[2]]["TP",]/(true_conf_matrix[[2]]["TP",]+true_conf_matrix[[2]]["FN",])
  # vioplot(Recall_L, Recall_V,col=rgb(0.1,0.4,0.7,0.7) , names=c("LASSO","VBSR"))
  # title(main="Recall")
  
  FPR_L<-confussion_matrices[[1]]["FP",]/(confussion_matrices[[1]]["FP",]+confussion_matrices[[1]]["TN",])
  FPR_V<-confussion_matrices[[2]]["FP",]/(confussion_matrices[[2]]["FP",]+confussion_matrices[[2]]["TN",])
  FPR_L_cons<-confussion_matrices_consensus[[1]]["FP",]/(confussion_matrices_consensus[[1]]["FP",]+confussion_matrices_consensus[[1]]["TN",])
  FPR_V_cons<-confussion_matrices_consensus[[2]]["FP",]/(confussion_matrices_consensus[[2]]["FP",]+confussion_matrices_consensus[[2]]["TN",])
  FPR_L_trueClust<-true_conf_matrix[[1]]["FP",]/(true_conf_matrix[[1]]["FP",]+true_conf_matrix[[1]]["TN",])
  FPR_V_trueClust<-true_conf_matrix[[2]]["FP",]/(true_conf_matrix[[2]]["FP",]+true_conf_matrix[[2]]["TN",])
  
  plot(FPR_L, Recall_L, type='p', xlim=c(0.0,0.2), ylim=c(0.5,1.0), xlab='FPR', ylab='TPR', pch=19,
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2);
  points(FPR_V, Recall_V, type='p', pch=17)
  #points(FPR_L_cons, Recall_L_cons, type='p', pch=19, col="red")
  #points(FPR_V_cons, Recall_V_cons, type='p', pch=17, col="red")
  points(FPR_L_trueClust, Recall_L_trueClust, type='p', pch=19, col="red")
  points(FPR_V_trueClust, Recall_V_trueClust, type='p', pch=17, col="red")
  title(main="TPR vs FPR")
  #legend("bottomright", 1.9, 
  #       c("LASSO","VBSR","All", "Consensus", "True Clusters" ), 
  #       col = c("black", "black", "black", "red", "green"),
  #       text.col = "black", lty = c(-1,-1, -1, -1, -1), pch = c(1, 2, 15, 15, 15),
  #       bg = "gray90")
  legend("bottomright", 1.9, 
         c("LASSO","VBSR","All", "True Clusters" ), 
         col = c("black", "black", "black", "red"),
         text.col = "black", lty = c(-1,-1, -1, -1), pch = c(1, 2, 15, 15),
         bg = "gray90",cex = 2)
  
  plot(Prec_L, Recall_L, type='p', xlim=c(0.5,1.0), ylim=c(0.5,1.0), xlab='Precision', ylab='Recall', pch=19,
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(Prec_V, Recall_V, type='p', pch=17)
  #points(Prec_L_cons, Recall_L_cons, type='p', pch=19, col="red")
  #points(Prec_V_cons, Recall_V_cons, type='p', pch=17, col="red")
  points(Prec_L_trueClust, Recall_L_trueClust, type='p', pch=19, col="red")
  points(Prec_V_trueClust, Recall_V_trueClust, type='p', pch=17, col="red")
  title(main="Precision vs Recall")
  # legend("bottomright", 1.9, 
  #        c("LASSO","VBSR","All", "Consensus", "True Clusters" ), 
  #        col = c("black", "black", "black", "red", "green"),
  #        text.col = "black", lty = c(-1,-1, -1, -1, -1), pch = c(1, 2, 15, 15, 15),
  #        bg = "gray90")
  legend("bottomleft", 1.9, 
         c("LASSO","VBSR","All", "True Clusters" ), 
         col = c("black", "black", "black", "red"),
         text.col = "black", lty = c(-1,-1, -1, -1), pch = c(1, 2, 15, 15),
         bg = "gray90",cex = 2)
  
  
  h<-hist(unlist(RsquareAjusted[[1]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  plot(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0.7,1), lty = "solid", xlab="Adjusted Rsquared",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  h<-hist(unlist(RsquareAjusted[[2]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
  
  # h<-hist(unlist(RsquareAjusted_consensus[[1]]), plot = FALSE, breaks = "FD")
  # cdf<-cumsum(h$counts)/sum(h$counts)
  # lines(h$mids, cdf, type = 'l', col = "blue2", lty = "dashed")
  # h<-hist(unlist(RsquareAjusted_consensus[[2]]), plot = FALSE, breaks = "FD")
  # cdf<-cumsum(h$counts)/sum(h$counts)
  # lines(h$mids, cdf, type = 'l', col = "red2", lty = "dashed")
  
  h<-hist(unlist(sim_results$AdjR2$LASSO), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  lines(h$mids, cdf, type = 'l', col = "blue", lty = "dashed")
  h<-hist(unlist(sim_results$AdjR2$VBSR), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  lines(h$mids, cdf, type = 'l', col = "red", lty = "dashed")
  
  # legend("topleft", 1.9, 
  #        c("LASSO","VBSR","All", "Consensus", "True Clusters" ), 
  #        col = c("blue", "red", "black", "black", "black"),
  #        text.col = "black", lty = c(1,1, 1, 2, 3),
  #        bg = "gray90")
  
  legend("topleft", 1.9, 
         c("LASSO","VBSR","All", "True Clusters" ), 
         col = c("blue", "red", "black", "black"),
         text.col = "black", lty = c(1,1, 1, 2),
         bg = "gray90", cex = 2)
  
  ##### Variance Explained ########
  h<-hist(unlist(VarEx[[1]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  plot(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,0.4), lty = "solid", xlab="Percentage of variance explained",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  h<-hist(unlist(VarEx[[2]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
  
  h<-hist(unlist(VarEx_consensus[[1]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  #lines(h$mids, cdf, type = 'l', col = "blue2", lty = "dashed")
  h<-hist(unlist(VarEx_consensus[[2]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  #lines(h$mids, cdf, type = 'l', col = "red2", lty = "dashed")
  
  h<-hist(unlist(sim_results$VarEx), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  lines(h$mids, cdf, type = 'l', col = "black", lty = "dashed")
  
  # legend("topleft", 1.9, 
  #        c("LASSO","VBSR","All", "Consensus", "True Clusters" ), 
  #        col = c("blue", "red", "black", "black", "black"),
  #        text.col = "black", lty = c(1,1, 1, 2, 3),
  #        bg = "gray90")
  
  legend("bottomright", 1.9, 
         c("LASSO","VBSR","All", "True Clusters" ), 
         col = c("blue", "red", "black", "black"),
         text.col = "black", lty = c(1,1, 1, 2),
         bg = "gray90", cex = 2)
  
  dev.off()
}

evaluate_sim_data<-function(sim, used_method="LINKER"){
  
  
  true_reg<-lapply(sim, function(x)x[[4]])
  NrSims<-length(sim)
  computed_betas_lasso<-list()
  computed_betas_VBSR<-list()
  confussion_matrix<-list()
  confussion_matrix[[1]]<-matrix(nrow = 4, ncol = NrSims)
  rownames(confussion_matrix[[1]])<-c("TP", "FP", "TN", "FN")
  confussion_matrix[[2]]<-matrix(nrow = 4, ncol = NrSims)
  rownames(confussion_matrix[[2]])<-c("TP", "FP", "TN", "FN")
  RsquareAjusted_lasso<-list()
  RsquareAjusted_vbsr<-list()
  VarEx<-list()
  
  for(i in 1:NrSims){
    Data<-sim[[i]][[1]][sim[[i]][[2]] , ]
    RegulatorData<-sim[[i]][[1]][sim[[i]][[3]] , ]
    Clusters<-sim[[i]][[5]]
    NrClusters<-length(unique(Clusters))
    computed_betas_lasso[[i]]<-matrix(ncol = NrClusters, nrow = nrow(RegulatorData))
    computed_betas_VBSR[[i]]<-matrix(ncol = NrClusters, nrow = nrow(RegulatorData))
    RsquareAjusted_lasso[[i]]<-numeric(length = NrClusters)
    RsquareAjusted_vbsr[[i]]<-numeric(length = NrClusters)
    VarEx[[i]]<-numeric(length = NrClusters)
    for (clustIdx in 1:NrClusters){
      
      CurrentClusterPositions = which(Clusters==clustIdx)
      nrGenesInClusters = length(CurrentClusterPositions)
      gene_svd<-svd(Data[CurrentClusterPositions,])
      
      if(used_method=="LINKER"){
        
        y<-gene_svd$v[,1]
        if(cor(apply(Data[CurrentClusterPositions,],2,mean), y)<0){
          y<- -y
        }
        scale(y, center = FALSE)
      }
      else{
        y<-apply(Data[CurrentClusterPositions,],2,mean)
      }
      
      X = RegulatorData
      
      fit = cv.glmnet(t(X), y, alpha = 1-1e-06, pmax=10)
      #fit = cv.glmnet(t(X), y, alpha = alpha, pmax = pmax, lambda = grid)
      
      nonZeroLambdas <- fit$lambda[which(fit$nzero>0)]
      nonZeroCVMs <- fit$cvm[which(fit$nzero>0)]
      
      if(length(which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE)))==0){
        
        #for now: just print a warning, *although* this error WILL cause LINKER to crash in a few steps.
        warnMessage <- paste0("\nOn cluster ",clustIdx," there were no cv.glm results that gave non-zero coefficients.")
        warning(warnMessage);
        bestNonZeroLambda<-fit$lambda.min
      }
      else{
        bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE))]
      }
      #b_o = coef(fit,s = bestNonZeroLambda)
      b_o = coef(fit,s = fit$lambda.min)
      b_opt_lasso <- c(b_o[2:length(b_o)]) # removing the intercept.
      
      res<-vbsr(y,t(X),n_orderings = 15,family='normal')
      betas<-res$beta
      max_beta<-max(abs(betas))
      #betas[which(abs(betas)<0.1*max_beta)]<-0
      betas[res$pval > 0.05/(nrow(RegulatorData)*nrow(Data))]<-0
      b_opt_vbsr<-betas
      
      VarEx[[i]][clustIdx]<-gene_svd$d[1]^2/sum(gene_svd$d^2)
      RsquareAjusted_lasso[[i]][clustIdx]<-compute_adj_r2(betas=b_opt_lasso, y, regs=X)
      RsquareAjusted_vbsr[[i]][clustIdx]<-compute_adj_r2(betas=b_opt_vbsr, y, regs=X)
      
      computed_betas_lasso[[i]][,clustIdx]<-b_opt_lasso
      computed_betas_VBSR[[i]][,clustIdx]<-b_opt_vbsr
    }
    
    confussion_matrix[[1]]["TN",i]=length(which( (computed_betas_lasso[[i]] == 0) & (true_reg[[i]]==0) ))
    confussion_matrix[[2]]["TN",i]=length(which( (computed_betas_VBSR[[i]] == 0) & (true_reg[[i]]==0) ))
    
    confussion_matrix[[1]]["TP",i]=length(which( (computed_betas_lasso[[i]] != 0) & (true_reg[[i]] != 0) ))
    confussion_matrix[[2]]["TP",i]=length(which( (computed_betas_VBSR[[i]] != 0) & (true_reg[[i]] != 0) ))
    
    confussion_matrix[[1]]["FP",i]=length(which( (computed_betas_lasso[[i]] != 0) & (true_reg[[i]]==0) ))
    confussion_matrix[[2]]["FP",i]=length(which( (computed_betas_VBSR[[i]] != 0) & (true_reg[[i]]==0) ))
    
    confussion_matrix[[1]]["FN",i]=length(which( (computed_betas_lasso[[i]] == 0) & (true_reg[[i]]!=0) ))
    confussion_matrix[[2]]["FN",i]=length(which( (computed_betas_VBSR[[i]] == 0) & (true_reg[[i]]!=0) ))
    
  }
  
  RsquareAjusted<-list(LASSO=RsquareAjusted_lasso, VBSR=RsquareAjusted_vbsr)
  
  sim_lasso<-unlist(computed_betas_lasso)
  sim_vbsr<-unlist(computed_betas_VBSR)
  true_reg<-unlist(true_reg)
  
  #scale
  true_betas<-scale(true_reg, center = FALSE)
  sim_betas_lasso<-scale(sim_lasso, center = FALSE)
  sim_betas_vbsr<-scale(sim_vbsr, center = FALSE)
  #
  
  TN_l=which( (sim_betas_lasso == 0) & (true_betas==0) )
  TN_v=which( (sim_betas_vbsr == 0) & (true_betas==0) )
  
  TP_l=which( (sim_betas_lasso != 0) & (true_betas != 0) )
  TP_v=which( (sim_betas_vbsr != 0) & (true_betas != 0) )
  
  FP_l=which( (sim_betas_lasso != 0) & (true_betas==0) )
  FP_v=which( (sim_betas_vbsr != 0) & (true_betas==0) )
  
  FN_l=which( (sim_betas_lasso == 0) & (true_betas!=0) )
  FN_v=which( (sim_betas_vbsr == 0) & (true_betas!=0) )
  
  model_all<-lm(true_betas[TP_l] ~ sim_betas_lasso[TP_l])
  plot(sim_betas_lasso[TP_l],true_betas[TP_l], pch=20, col="blue",
       xlab= expression(paste("true weights ", beta)), ylab=expression(paste("estimated weights ", beta)),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2);
  abline(model_all$coefficients, col="darkblue");
  points(sim_betas_lasso[FP_l],true_betas[FP_l], pch=17, col="purple")
  points(sim_betas_lasso[FN_l],true_betas[FN_l], pch=16, col="green")
  legend("topleft", 1.9, 
         c(paste0("MSE: ", round(sum((true_betas-sim_betas_lasso)^2)/length(true_betas),digits = 4)),paste0("TP: ",length(TP_l)), paste0("FP: ",length(FP_l)), paste0("FN: ",length(FN_l))), 
         col = c(NA,"blue", "purple", "green"),
         text.col = "black", lty = c(-1,-1, -1, -1), pch = c(NA,20, 17, 16),
         merge = TRUE, bg = "gray90", cex = 2)
  abline(a=0,b=1, col="black")
  
  model_all<-lm(true_betas[TP_v] ~ sim_betas_vbsr[TP_v])
  plot(sim_betas_vbsr[TP_v],true_betas[TP_v], pch=20, col="blue",
       xlab= expression(paste("true weights ", beta)), ylab=expression(paste("estimated weights ", beta)),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2);
  abline(model_all$coefficients, col="darkblue");
  points(sim_betas_vbsr[FP_v],true_betas[FP_v], pch=17, col="purple")
  points(sim_betas_vbsr[FN_v],true_betas[FN_v], pch=16, col="green")
  legend("topleft", 1.9, 
          c(paste0("MSE: ",round(sum((true_betas-sim_betas_vbsr)^2)/length(true_betas),digits = 4)), paste0("TP: ",length(TP_v)), paste0("FP: ",length(FP_v)), paste0("FN: ",length(FN_v))), 
          col = c(NA,"blue", "purple", "green"),
          text.col = "black", lty = c(-1,-1, -1, -1), pch = c(NA,20, 17, 16),
          merge = TRUE, bg = "gray90", cex = 2)
  abline(a=0,b=1, col="black")

  
  return(list(conf_matrix=confussion_matrix,AdjR2=RsquareAjusted, VarEx=VarEx))
}

compute_adj_r2<-function(betas, y, regs)
{
  prediction<-betas %*% regs
  NrRegulators<-length(which(betas != 0))
  SStot = sum((y-mean(y))^2)
  SSres = sum( (prediction-y)^2)
  Rsquare = 1 - (SSres / SStot)
  return( 1 - ( ( (ncol(regs) - 1) / (ncol(regs) - NrRegulators) ) * (1-Rsquare) ) )
}

plot_res_real_data<-function(res){
  
  pdf("plots_real_data.pdf")
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0)) 
  
  NrSims<-length(res[[1]])
  
  
  Nr_bootstraps<-length(res[[1]][[1]]$bootstrapResults)

  RsquareAjusted<-list()
  RsquareAjusted_test<-list()
  RsquareAjusted_consensus<-list()
  VarEx<-list()
  VarEx_test<-list()
  VarEx_consensus<-list()
  
  
  for(idx in 1:2){
    all_idx<-1
    RsquareAjusted[[idx]]<-list(length = NrSims*Nr_bootstraps)
    RsquareAjusted_test[[idx]]<-list(length = NrSims*Nr_bootstraps)
    RsquareAjusted_consensus[[idx]]<-list(length = NrSims)
    VarEx[[idx]]<-list(length = NrSims*Nr_bootstraps)
    VarEx_test[[idx]]<-list(length = NrSims*Nr_bootstraps)
    VarEx_consensus[[idx]]<-list(length = NrSims)
    
    for(j in 1:NrSims){
      
      for(i in 1: Nr_bootstraps){
        
        RsquareAjusted[[idx]][[all_idx]]<-rep(res[[idx]][[j]]$bootstrapResults[[i]]$trainingStats[,"RsquareAdjusted"],res[[idx]][[j]]$bootstrapResults[[i]]$trainingStats[,"nrGen"])
        VarEx[[idx]][[all_idx]]<-rep(res[[idx]][[j]]$bootstrapResults[[i]]$trainingStats[,"condition"],res[[idx]][[j]]$bootstrapResults[[i]]$trainingStats[,"nrGen"])
        RsquareAjusted_test[[idx]][[all_idx]]<-rep(res[[idx]][[j]]$bootstrapTestStats[[i]][,"RsquareAdjusted"],res[[idx]][[j]]$bootstrapResults[[i]]$trainingStats[,"nrGen"])
        VarEx_test[[idx]][[all_idx]]<-rep(res[[idx]][[j]]$bootstrapTestStats[[i]][,"condition"],res[[idx]][[j]]$bootstrapResults[[i]]$trainingStats[,"nrGen"])
        
        # Clusters<-res[[idx]][[j]]$bootstrapResults[[i]]$ModuleMembership
        # NrClusters<-res[[idx]][[j]]$bootstrapResults[[i]]$NrModules
        # 
        # GeneNames<-res[[idx]][[j]]$bootstrapResults[[i]]$AllGenes
        # 
        # silhouette<-numeric(length = length(GeneNames))
        # for(geneIdx in 1:length(GeneNames)){
        #   
        #   geneName<-GeneNames[geneIdx]
        #   geneCluster<-Clusters[geneName,]
        #   
        #   if(geneName %in% rownames(Data) == FALSE){
        #     print("foo")
        #   }
        #   gene<-Data[geneName,]
        #   
        #   max_outClust_cor<-0
        #   for (clustIdx in 1:NrClusters){
        #     
        #     CurrentClusterPositions = which(Clusters==clustIdx)
        #     nrGenesInClusters = length(CurrentClusterPositions)
        #     module_data<-Data[rownames(Clusters)[CurrentClusterPositions],]
        #     avg_corr<-(sum(abs(cor(gene,t(module_data)))) - 1)/nrow(module_data) # We need to remove the corr with itself
        #     if(clustIdx==geneCluster){
        #       inClust_dis<-1-avg_corr
        #     }
        #     else if(avg_corr > max_outClust_cor){
        #       outClust_dis<-1-avg_corr
        #       max_outClust_cor<-avg_corr
        #     }
        #   }
        #   silhouette[geneIdx]<-(outClust_dis-inClust_dis)/max(outClust_dis,inClust_dis)
        #   
        # }
       
        
        all_idx<-all_idx+1
      }
      
      #RsquareAjusted_consensus[[idx]][[j]]<-rep(res[[idx]][[j]]$consensusTestStats[,"RsquareAdjusted"],res[[idx]][[j]]$consensusTestStats[,"nrGen"])
      #VarEx_consensus[[idx]][[j]]<-rep(res[[idx]][[j]]$consensusTestStats[,"condition"],res[[idx]][[j]]$consensusTestStats[,"nrGen"])
    }
  }
  
  
  h<-hist(unlist(RsquareAjusted[[1]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  plot(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "solid", xlab="Adjusted Rsquared",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  h<-hist(unlist(RsquareAjusted[[2]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
  
  # h<-hist(unlist(RsquareAjusted_consensus[[1]]), plot = FALSE, breaks = "FD")
  # cdf<-cumsum(h$counts)/sum(h$counts)
  # lines(h$mids, cdf, type = 'l', col = "blue2", lty = "dashed")
  # h<-hist(unlist(RsquareAjusted_consensus[[2]]), plot = FALSE, breaks = "FD")
  # cdf<-cumsum(h$counts)/sum(h$counts)
  # lines(h$mids, cdf, type = 'l', col = "red2", lty = "dashed")
  legend("topleft", 1.9, 
         c("LASSO","VBSR","Training", "Test"  ), 
         col = c("blue", "red", "black", "black"),
         text.col = "black", lty = c(1,1, 1, 2),
         bg = "gray90", cex = 2)
  
  all_idx<-1
  for(j in 1:NrSims){
    
    for(i in 1: Nr_bootstraps){
      h<-hist(RsquareAjusted[[1]][[all_idx]], plot = FALSE, breaks = "FD")
      cdf<-cumsum(h$counts)/sum(h$counts)
      lines(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "solid", xlab="Adjusted Rsquared")
      h<-hist(RsquareAjusted[[2]][[all_idx]], plot = FALSE, breaks = "FD")
      cdf<-cumsum(h$counts)/sum(h$counts)
      lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
      
      h<-hist(RsquareAjusted_test[[1]][[all_idx]], plot = FALSE, breaks = "FD")
      cdf<-cumsum(h$counts)/sum(h$counts)
      lines(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "dotted", xlab="Adjusted Rsquared")
      h<-hist(RsquareAjusted_test[[2]][[all_idx]], plot = FALSE, breaks = "FD")
      cdf<-cumsum(h$counts)/sum(h$counts)
      lines(h$mids, cdf, type = 'l', col = "red4", lty = "dotted")
      all_idx<-all_idx+1
    }
    # h<-hist(RsquareAjusted_consensus[[1]][[j]], plot = FALSE, breaks = "FD")
    # cdf<-cumsum(h$counts)/sum(h$counts)
    # lines(h$mids, cdf, type = 'l', col = "blue2", lty = "dashed")
    # h<-hist(RsquareAjusted_consensus[[2]][[j]], plot = FALSE, breaks = "FD")
    # cdf<-cumsum(h$counts)/sum(h$counts)
    # lines(h$mids, cdf, type = 'l', col = "red2", lty = "dashed")
  }
  

  ##### Variance Explained ########
  h<-hist(unlist(VarEx[[1]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  plot(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "solid", xlab="Percentage of variance explained",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  h<-hist(unlist(VarEx[[2]]), plot = FALSE, breaks = "FD")
  cdf<-cumsum(h$counts)/sum(h$counts)
  lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
  
  # h<-hist(unlist(VarEx_consensus[[1]]), plot = FALSE, breaks = "FD")
  # cdf<-cumsum(h$counts)/sum(h$counts)
  # lines(h$mids, cdf, type = 'l', col = "blue2", lty = "dashed")
  # h<-hist(unlist(VarEx_consensus[[2]]), plot = FALSE, breaks = "FD")
  # cdf<-cumsum(h$counts)/sum(h$counts)
  # lines(h$mids, cdf, type = 'l', col = "red2", lty = "dashed")
  
  
  legend("bottomright", 1.9, 
         c("LASSO","VBSR","Training", "Test" ), 
         col = c("blue", "red", "black", "black"),
         text.col = "black", lty = c(1,1, 1, 2),
         bg = "gray90", cex = 2)
  
  all_idx<-1
  for(j in 1:NrSims){
    
    for(i in 1: Nr_bootstraps){
      h<-hist(VarEx[[1]][[all_idx]], plot = FALSE, breaks = "FD")
      cdf<-cumsum(h$counts)/sum(h$counts)
      lines(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "dotted", xlab="Adjusted Rsquared")
      h<-hist(VarEx[[2]][[all_idx]], plot = FALSE, breaks = "FD")
      cdf<-cumsum(h$counts)/sum(h$counts)
      lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
      
      h<-hist(VarEx_test[[1]][[all_idx]], plot = FALSE, breaks = "FD")
      cdf<-cumsum(h$counts)/sum(h$counts)
      lines(h$mids, cdf, type = 'l', col = "blue4", xlim = c(0,1), lty = "dotted", xlab="Adjusted Rsquared")
      h<-hist(VarEx_test[[2]][[all_idx]], plot = FALSE, breaks = "FD")
      cdf<-cumsum(h$counts)/sum(h$counts)
      lines(h$mids, cdf, type = 'l', col = "red4", lty = "solid")
      all_idx<-all_idx+1
    }
    # h<-hist(VarEx_consensus[[1]][[j]], plot = FALSE, breaks = "FD")
    # cdf<-cumsum(h$counts)/sum(h$counts)
    # lines(h$mids, cdf, type = 'l', col = "blue2", lty = "dashed")
    # h<-hist(VarEx_consensus[[2]][[j]], plot = FALSE, breaks = "FD")
    # cdf<-cumsum(h$counts)/sum(h$counts)
    # lines(h$mids, cdf, type = 'l', col = "red2", lty = "dashed")
  }
  dev.off()
}

