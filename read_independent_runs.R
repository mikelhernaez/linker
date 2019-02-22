#OV
link_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM")
GEAs<-list()
graphs<-list()
for(i in 1:length(link_mode)){
  tmp<-readRDS(paste0("/data/blatti/linker_TCGA/linker_190210/linker-Tumor_OV50-",link_mode[i],"/Tumor_OV50.tar8855_reg638.m100_b10.rds"))
  graph_mode<-names(tmp$graphs[[link_mode[i]]])
  GEAs[[ link_mode[i] ]]<-list()
  graphs[[ link_mode[i] ]]<-list()
  for(j in 1:length(graph_mode)){
    GEAs[[ link_mode[i] ]][[ graph_mode[j] ]] <- tmp$GEAs[[ link_mode[i] ]][[ graph_mode[j] ]]
    graphs[[ link_mode[i] ]][[ graph_mode[j] ]] <- tmp$graphs[[ link_mode[i] ]][[ graph_mode[j] ]]
    print(paste0("Done for (",link_mode[i],",",graph_mode[j], ") computed!"))    
  }
}

graph_mode<-c("VBSR", "LASSOmin", "LASSO1se", "LM")
GEAs$NET<-list()
graphs$NET<-list()
for(j in 1:length(graph_mode)){
  tmp<-readRDS(paste0("/data/blatti/linker_TCGA/linker_190210/linker-Tumor_OV50-",graph_mode[j],"/Tumor_OV50.tar8855_reg638.single_gene.rds"))
  GEAs$NET[[ graph_mode[j] ]] <- tmp$GEAs[[ graph_mode[j] ]]
  graphs$NET[[ graph_mode[j] ]] <- tmp$graphs[[ graph_mode[j] ]]
  print(paste0("Done for (NET,",graph_mode[j], ") computed!"))    
}

AUC_regs_edge<-LINKER_plot_GEAs_norm_regs(GEAs, graphs,type="EDGE", max_y = c(0.3,0.6,0.6,2), min_x = c(-20,-50,-40,-50))
AUC_regs_reg<-LINKER_plot_GEAs_norm_regs(GEAs, graphs,type="REG", max_y = c(.01,.03,.05,.1), min_x = c(-20,-50,-40,-50))

AUC_sets_edge<-LINKER_plot_GEAs_norm_sets(GEAs, graphs,type="EDGE", max_y = c(0.3,2,0.6,20), min_x = c(-20,-50,-40,-50))
AUC_sets_reg<-LINKER_plot_GEAs_norm_sets(GEAs, graphs,type="REG", max_y = c(0.01,.15,0.02,0.5), min_x = c(-20,-50,-40,-50))

AUC_regsets_edge<-LINKER_plot_GEAs_norm_regsets(GEAs, graphs,type="EDGE", max_y = c(2,15,4,60), min_x = c(-20,-50,-40,-50))
AUC_regsets_reg<-LINKER_plot_GEAs_norm_regsets(GEAs, graphs,type="REG", max_y = c(.1,2,.2,4), min_x = c(-20,-50,-40,-50))

ggplot(AUC_regs_edge, aes(link_mode, Rank)) + geom_boxplot(aes(fill=DB)) + scale_y_reverse() +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),width = .9, linetype = "dashed", color="red")+
  #theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="Rank of GEAreg elements per 1K graph edges", 
       x="Method for generating the modules",
       y="Rank")

ggplot(AUC_regs_reg, aes(link_mode, Rank)) + geom_boxplot(aes(fill=DB)) + 
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),width = .9, linetype = "dashed", color="red")+
  labs(title="Rank of GEAreg elements per graph regulator", 
       x="Method for generating the modules",
       y="Rank")

ggplot(AUC_sets_edge, aes(link_mode, Rank)) + geom_boxplot(aes(fill=DB)) + 
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),width = .9, linetype = "dashed", color="red")+
  labs(title="Rank of GEAset elements per 1K graph edges", 
       x="Method for generating the modules",
       y="Rank")

ggplot(AUC_sets_reg, aes(link_mode, Rank)) + geom_boxplot(aes(fill=DB)) + 
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),width = .9, linetype = "dashed", color="red")+
  labs(title="Rank of GEAset elements per graph regulator", 
       x="Method for generating the modules",
       y="Rank")

ggplot(AUC_regsets_edge, aes(link_mode, Rank)) + geom_boxplot(aes(fill=DB)) + 
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),width = .9, linetype = "dashed", color="red")

ggplot(AUC_regsets_reg, aes(link_mode, Rank)) + geom_boxplot(aes(fill=DB)) + 
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),width = .9, linetype = "dashed", color="red")



#geom_point(data=AUC_regs_edge %>% group_by(link_mode) %>% summarise(ave = mean(Rank)),aes(link_mode,ave)) + 
#HNSC
link_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM")
GEAs<-list()
graphs<-list()
for(i in 1:length(link_mode)){
  tmp<-readRDS(paste0("/data/blatti/linker_TCGA/linker_190210/linker-Tumor_HNSC50-",link_mode[i],"/Tumor_HNSC50.tar8791_reg702.m100_b10.rds"))
  graph_mode<-names(tmp$graphs[[link_mode[i]]])
  GEAs[[ link_mode[i] ]]<-list()
  graphs[[ link_mode[i] ]]<-list()
  for(j in 1:length(graph_mode)){
    GEAs[[ link_mode[i] ]][[ graph_mode[j] ]] <- tmp$GEAs[[ link_mode[i] ]][[ graph_mode[j] ]]
    graphs[[ link_mode[i] ]][[ graph_mode[j] ]] <- tmp$graphs[[ link_mode[i] ]][[ graph_mode[j] ]]
    print(paste0("Done for (",link_mode[i],",",graph_mode[j], ") computed!"))    
  }
}
  
graph_mode<-c("VBSR", "LASSOmin", "LASSO1se", "LM")
GEAs$NET<-list()
graphs$NET<-list()
for(j in 1:length(graph_mode)){
  tmp<-readRDS(paste0("/data/blatti/linker_TCGA/linker_190210/linker-Tumor_HNSC50-",graph_mode[j],"/Tumor_HNSC50.tar8791_reg702.single_gene.rds"))
  GEAs$NET[[ graph_mode[j] ]] <- tmp$GEAs[[ graph_mode[j] ]]
  graphs$NET[[ graph_mode[j] ]] <- tmp$graphs[[ graph_mode[j] ]]
  print(paste0("Done for (NET,",graph_mode[j], ") computed!"))    
}




LINKER_plot_GEAs_norm_regs(GEAs, graphs,type="EDGE", max_y = c(0.3,0.6,0.6,2), min_x = c(-20,-50,-40,-50))
LINKER_plot_GEAs_norm_regs(GEAs, graphs,type="REG", max_y = c(.01,.03,.05,.1), min_x = c(-20,-50,-40,-50))

LINKER_plot_GEAs_norm_sets(GEAs, graphs,type="EDGE", max_y = c(0.3,0.6,0.6,2), min_x = c(-20,-50,-40,-50))
LINKER_plot_GEAs_norm_sets(GEAs, graphs,type="REG", max_y = c(.01,.03,.05,.1), min_x = c(-20,-50,-40,-50))

LINKER_plot_GEAs_norm_regsets(GEAs, graphs,type="EDGE", max_y = c(0.3,0.6,0.6,2), min_x = c(-20,-50,-40,-50))
LINKER_plot_GEAs_norm_regsets(GEAs, graphs,type="REG", max_y = c(.1,.3,.5,1), min_x = c(-20,-50,-40,-50))
