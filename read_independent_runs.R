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

draw_all(GEAs, graphs, FDR=0.05)

# 
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

draw_all(GEAs, graphs, FDR=0.05)
