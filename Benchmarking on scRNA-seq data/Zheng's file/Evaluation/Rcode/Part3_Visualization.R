#part4 visulization
rm(list = ls())

##1. Loading packages
##--------------------------------------------------------------------------------------------
library(ggplot2)
library(RColorBrewer)
##--------------------------------------------------------------------------------------------

#2. Set up the overall file path
##--------------------------------------------------------------------------------------------
data_source<-"Zheng"
dim_data<-300
p_TF<-8
#2. Set up the overall file path
##--------------------------------------------------------------------------------------------
File_path<-paste("./Benchmarking on scRNA-seq data/",data_source,"'s file",sep = "")

##--------------------------------------------------------------------------------------------
##3. Load data
##--------------------------------------------------------------------------------------------
load(file = paste(File_path,"/Evaluation/Obj_use/AUPR_ratio_res_densy_sorted_right.Rdata",sep = ""))
##--------------------------------------------------------------------------------------------
AUPR_ratio_res_densy<-AUPR_ratio_res_densy_sorted
##4. Set the order of celltype and method
##--------------------------------------------------------------------------------------------
celltype_reorder<-c("B cells","Monocytes","CD8+ T cells","CD4+ T cells")
method_order<-c(row.names(AUPR_ratio_res_densy))
save(method_order,file = paste(File_path,"/Evaluation/Obj_use/method_order.Rdata",sep = ""))
##--------------------------------------------------------------------------------------------

##5. Description of data
##--------------------------------------------------------------------------------------------
load(paste(File_path,"/Evaluation/Obj_use/description_data_4group.Rdata",sep = ""))
celltypes_num<-as.numeric(description_data["celltypes_num"])
p_TF<-as.numeric(description_data["p_TF"])
p_nonTF<-as.numeric(description_data["p_nonTF"])
p_GRN<-as.numeric(description_data["p_GRN"])
celltype_order<-colnames(AUPR_ratio_res_densy)

##6. Visualization for each criterion
##--------------------------------------------------------------------------------------------
##6.1 AUPR
AUPR_ratio_res_densy_nor<-AUPR_ratio_res_densy
AUPR_ratio_res_densy_nor<-sapply(X = 1:celltypes_num,FUN = function(g){return((AUPR_ratio_res_densy[,g])/(max(AUPR_ratio_res_densy[,g])))})
colnames(AUPR_ratio_res_densy_nor)<-colnames(AUPR_ratio_res_densy)
AUPR_ratio_res_densy_nor<-ifelse(AUPR_ratio_res_densy >=1,AUPR_ratio_res_densy_nor,NA)
##
plotdata_AUPR_densy<-data.frame(Method = rep(rownames(AUPR_ratio_res_densy_nor),ncol(AUPR_ratio_res_densy_nor)),
                     Celltype = rep(colnames(AUPR_ratio_res_densy_nor),each = nrow(AUPR_ratio_res_densy_nor)),
                     AUPR_ratio = as.vector(AUPR_ratio_res_densy_nor))
plotdata_AUPR_densy$Method<-factor(plotdata_AUPR_densy$Method,levels = rev(method_order))
plotdata_AUPR_densy$Celltype<-factor(plotdata_AUPR_densy$Celltype,levels = celltype_reorder)
plotdata_AUPR_densy$AUPR_ratio<-as.numeric(plotdata_AUPR_densy$AUPR_ratio)

anno_data<-matrix(NA,nrow = 0,ncol = 3)
for(g in 1:celltypes_num){
  num_method_nonNA<-sum(!is.na(AUPR_ratio_res_densy_nor[,g]))
  anno_data<-rbind(anno_data,  cbind(rep(celltype_order[g],num_method_nonNA),
                                     names(AUPR_ratio_res_densy[!is.na(AUPR_ratio_res_densy_nor[,g]),g]),
                                     as.vector(round(AUPR_ratio_res_densy[!is.na(AUPR_ratio_res_densy_nor[,g]),g],1))))
}
anno_data<-as.data.frame(anno_data)
names(anno_data)<-c("Celltype","Method","AUPR_ratio")
anno_data$Celltype<-factor(anno_data$Celltype,levels = celltype_reorder)
anno_data$Method<-factor(anno_data$Method,levels = rev(method_order))
anno_data$AUPR_ratio<-as.numeric(anno_data$AUPR_ratio)

p_AUPR_densy <- ggplot(data = plotdata_AUPR_densy, aes(x = Celltype, y = Method)) +
  geom_tile(aes(fill = AUPR_ratio), width = 0.95) +
  scale_fill_gradientn(colours = c("#fb8353", "#ffcfb8","#fbf3ed"),
                       values = c(1, 0.5, 0),
                       na.value = "#C0BBBB",
                       limits = c(0, 1),
                       breaks = c(0, 1),
                       labels = c("Low", "High"),
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  geom_text(data = anno_data, aes(x = Celltype, y = Method, label = AUPR_ratio), size = 4.1, cex = 6*4) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 0, face = "bold"),  # Bold x-axis labels
        axis.text.y = element_text(face = "bold"), 
        legend.position = "right",
        legend.key.height = unit(20, "pt"),
        legend.key.width = unit(20, "pt"),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.border = element_blank()) +
  guides(color = FALSE) +
  labs(fill = "AUPRC ratio")

p_AUPR_densy

ggsave(paste(File_path,"/Evaluation/picture_heatmap/AUPR_Zheng.pdf",sep = ""),
       width = 5.5, height = 4, scale = 1, dpi = 300, bg = "transparent")
