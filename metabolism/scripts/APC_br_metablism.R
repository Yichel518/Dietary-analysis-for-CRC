setwd("F:/A_MetaboAnalysis/Black_rice/Black_rice_metab/")
#################################### Apc BR #########################################
library(readxl)
library(tidyverse)
library(limma)
ALL_sample_data <- read_excel("F:/A_MetaboAnalysis/Black_rice/Black_rice_metab/Apc_results/1.Data_Assess/all_group/ALL_sample_data.xlsx")
colnames(ALL_sample_data) <- gsub("_1_","_14w_",colnames(ALL_sample_data))%>%
  gsub("_2_","_22w_",.)
colnames(ALL_sample_data)
QC <- ALL_sample_data[,168:184]
#View(ALL_sample_data)
#ALL_meta <- data.frame(Compounds=ALL_sample_data$Compounds,ALL_sample_data[,grep("_14w_",colnames(ALL_sample_data))])
ALL_meta <- data.frame(Compounds=ALL_sample_data$Compounds,ALL_sample_data[,12:167])
colnames(ALL_meta)
ALL_meta[1:5,1:5]
rownames(ALL_meta)<- ALL_meta$Compounds
ALL_meta$Compounds <- NULL
head(name)
colnames(ALL_meta)
ALL_meta <- ALL_meta[,-c(17:48,95:126)]
sampleInfo <- data.frame(strsplit2(colnames(ALL_meta),"_"))
colnames(sampleInfo) <- c("condition","sex","stage","replicate")
sampleInfo$group <- paste0(sampleInfo$condition,"_",sampleInfo$sex)
sampleInfo$group1 <- paste0(sampleInfo$condition,"_",sampleInfo$stage)
sampleInfo$group2 <- paste0(sampleInfo$condition,"_",sampleInfo$sex,"_",sampleInfo$stage)
rownames(sampleInfo) <- colnames(ALL_meta)
sampleInfo$group1 <- factor(sampleInfo$group1,levels = c("WT_14w","C_14w","B_14w",
                                                         "WT_22w","C_22w","B_22w"))
dim(sampleInfo)
colnames(ALL_meta)
logc <- log(ALL_meta+1)
head(logc)
rownames(ALL_meta)[1:5]
dim(logc)
sampleInfo

scale_logc <- t(scale(t(logc)))
library(PCAtools)
pca <- pca(logc, metadata =sampleInfo)
biplot(pca, x = 'PC1', y = 'PC2') #可以看到这两个成分对样品的解释度
#screen(p) #可以看到所有成分对样品的解释
##将pca与关联样本信息表 
pca$variance
pca_rotated_plus <- cbind(pca$rotated,sampleInfo)
library(ggsci)
pca_rotated_plus$condition <- factor(pca_rotated_plus$condition,levels=c("WT","C","B"))
pca_rotated_plus$group <- factor(pca_rotated_plus$group1,levels = c("WT_14w","C_14w","B_14w","WT_22w","C_22w","B_22w"))
ggplot(pca_rotated_plus,aes(x = PC1 , y = PC2))+
  geom_point(size = 8,aes(#shape = group, 
    fill = condition,color = condition)) +
  stat_ellipse(aes(color = condition,fill=pca_rotated_plus$condition),linetype = 'dashed',size = 1, show.legend = FALSE) + #添加分组椭圆
  labs(x = 'PC1 (23.72%)',y = 'PC2 (14.89%)') + 
  3  scale_shape_manual(values = c(21,22,21,22,21,22,21,22,21,22))+ 
  scale_fill_manual(values = c("#1D5DC4","#F8C227","#0094A5"))+
  scale_color_manual(values = c("#1D5DC4","#F8C227","#0094A5"))+
  theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.position = "top",legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))
#"#2E8B57","#88ada6","#D8BFD8","#800000","#ffc773","#a1afc9","#4682B4","#808080","#556B2F","#5F9EA0"

ggplot(pca_rotated_plus,aes(x = PC1 , y = PC2))+
  geom_point(size = 8,aes(shape = group, fill = group)) +
  stat_ellipse(aes(color = group1,fill=group1),linetype = 'dashed',size = 1, show.legend = FALSE) + #添加分组椭圆
  labs(x = 'PC1 (34.29%)',y = 'PC2 (12.73%)') + 
  scale_shape_manual(values = c(21,21,21,22,22,22))+ 
  scale_fill_manual(values = c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))+
  theme(legend.title=element_blank())+theme(legend.position = 'top')

RSD <- data.frame(rsd=apply(QC,1,sd)/apply(QC,1,mean))
rownames(RSD) <- rownames(logc)
RSD0.2 <- subset(RSD,RSD$rsd< 0.2)
logc <- logc[rownames(RSD0.2),]
dim(logc)
dim(ALL_meta)
#View(logc)
ALL_meta <- ALL_meta[rownames(logc),]
dim(RSD)
colnames(logc)
scale_logc <- t(scale(t(logc)))


############################ pvalue FC C vs WT 22w  ################################
Pvalue<-c(rep(0,nrow(logc))) 
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._22w",colnames(logc))])==0&&sd(logc[i,grep("WT_._22w",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._22w",colnames(logc))]),as.numeric(logc[i,grep("WT_._22w",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._22w",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("WT_._22w",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)
head(mix)
mix$log2FC <- log2(as.numeric(mix$FC))
#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(63:92)])), 
                sampleInfo$condition[c(63:92)], orthoI = 1) 
vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
colnames(mix)
mix_select=subset(mix,FDR<0.05 & .>1 )
dim(mix_select)
colnames(mix)
res_MF<- mix
#View(res_MF)
res_MF$type <- ifelse(res_MF$FDR< 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')
res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
res_MF$type[which(res_MF$VIP<1)] <- "noSig"
mix_select <- mix_select[order(mix_select$log2FC),]
mix_select$group <- 0
mix_select$group[which(mix_select$log2FC>0)] <- "Control"
mix_select$group[which(mix_select$log2FC<0)] <- "WildType"
Diff_C_WT <- mix_select
#View(Diff_C_WT)
colnames(ALL_sample_data)
Diff_C_WT$Compounds <- rownames(Diff_C_WT)

######################## pvalue FC   C vs B 22w ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._22w",colnames(logc))])==0&&sd(logc[i,grep("B_._22w",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._22w",colnames(logc))]),as.numeric(logc[i,grep("B_._22w",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._22w",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("B_._22w",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(47:78)])), 
                sampleInfo$condition[c(47:78)], orthoI = 1) 
vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix,vip)
colnames(mix)
mix$log2FC <- log2(as.numeric(mix$FC))
mix_select=subset(mix,FDR<0.05 & .>1 )
dim(mix_select)
colnames(mix)
res_MF<- mix
#View(res_MF)
res_MF$type <- ifelse(res_MF$FDR< 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')

res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
res_MF$type[which(res_MF$VIP<1)] <- "noSig"
mix_select <- mix_select[order(mix_select$log2FC),]
mix_select$group <- 0
mix_select$group[which(mix_select$log2FC>0)] <- "Control"
mix_select$group[which(mix_select$log2FC<0)] <- "Blackrice"
Diff_C_B <- mix_select
Diff_C_B$Compounds <- rownames(Diff_C_B)
dim(Diff_C_B)


######################## pvalue FC   WT vs B 22w ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("WT_._22w",colnames(logc))])==0&&sd(logc[i,grep("B_._22w",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("WT_._22w",colnames(logc))]),as.numeric(logc[i,grep("B_._22w",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("WT_._22w",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("B_._22w",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(79:92,47:62)])), 
                sampleInfo$condition[c(79:92,47:62)], orthoI = 1) 
vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix,vip)
colnames(mix)
mix$log2FC <- log2(as.numeric(mix$FC))
mix_select=subset(mix,FDR<0.05 & .>1 )
dim(mix_select)
colnames(mix)
res_MF<- mix
#View(res_MF)
res_MF$type <- ifelse(res_MF$FDR< 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')

res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
res_MF$type[which(res_MF$VIP<1)] <- "noSig"
mix_select <- mix_select[order(mix_select$log2FC),]
mix_select$group <- 0
mix_select$group[which(mix_select$log2FC>0)] <- "WildType"
mix_select$group[which(mix_select$log2FC<0)] <- "Blackrice"
Diff_WT_B <- mix_select
Diff_WT_B$Compounds <- rownames(Diff_WT_B)
dim(Diff_WT_B)

bc=Diff_C_B[Diff_C_B$group=="Blackrice",]
cb=Diff_C_B[Diff_C_B$group=="Control",]
wtb=Diff_WT_B[Diff_WT_B$group=="WildType",]
bwt=Diff_WT_B[Diff_WT_B$group=="Blackrice",]
cwt=Diff_C_WT[Diff_C_WT$group=="Control",]
wtc=Diff_C_WT[Diff_C_WT$group=="WildType",]
h14 <- intersect(wtc$Compounds,bc$Compounds)
d14 <- intersect(cwt$Compounds,bwt$Compounds)


########################### pvalue FC C vs WT 14w  ################################
Pvalue<-c(rep(0,nrow(logc))) 
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._14w",colnames(logc))])==0&&sd(logc[i,grep("WT_._14w",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._14w",colnames(logc))]),as.numeric(logc[i,grep("WT_._14w",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._14w",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("WT_._14w",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)
head(mix)
mix$log2FC <- log2(as.numeric(mix$FC))
#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("C_._14w",colnames(scale_logc)),
                                              grep("WT_._14w",colnames(scale_logc)))])), 
                sampleInfo[c(grep("C_._14w",colnames(scale_logc)),
                             grep("WT_._14w",colnames(scale_logc))),]$condition,orthoI = 1) 
vip<- getVipVn(oplsda)%>%as.data.frame(.)
#View(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
colnames(mix)
#mix$log2FC <- log2(as.numeric(mix$FC))
mix_select=subset(mix,FDR<0.05 & .>1 )
dim(mix_select)
colnames(mix)
res_MF<- mix
#View(res_MF)
res_MF$type <- ifelse(res_MF$FDR< 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')
res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
res_MF$type[which(res_MF$VIP<1)] <- "noSig"
mix_select <- mix_select[order(mix_select$log2FC),]
mix_select$group <- 0

mix_select$group[which(mix_select$log2FC>0)] <- "Control"
mix_select$group[which(mix_select$log2FC<0)] <- "WildType"
Diff_C_WT <- mix_select
#View(Diff_C_WT)
colnames(ALL_sample_data)
Diff_C_WT$Compounds <- rownames(Diff_C_WT)

######################## pvalue FC   C vs B 14w ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._14w",colnames(logc))])==0&&sd(logc[i,grep("B_._14w",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._14w",colnames(logc))]),as.numeric(logc[i,grep("B_._14w",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._14w",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("B_._14w",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("C_._14w",colnames(scale_logc)),
                                              grep("B_._14w",colnames(scale_logc)))])), 
                sampleInfo[c(grep("C_._14w",colnames(scale_logc)),
                             grep("B_._14w",colnames(scale_logc))),]$condition,orthoI = 1) 
vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix,vip)
colnames(mix)
mix$log2FC <- log2(as.numeric(mix$FC))
mix_select=subset(mix,FDR<0.05 & .>1 )
dim(mix_select)
colnames(mix)
res_MF<- mix
#View(res_MF)
res_MF$type <- ifelse(res_MF$FDR< 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')

res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
res_MF$type[which(res_MF$VIP<1)] <- "noSig"
mix_select <- mix_select[order(mix_select$log2FC),]
mix_select$group <- 0
mix_select$group[which(mix_select$log2FC>0)] <- "Control"
mix_select$group[which(mix_select$log2FC<0)] <- "Blackrice"
Diff_C_B <- mix_select
Diff_C_B$Compounds <- rownames(Diff_C_B)
dim(Diff_C_B)

######################## pvalue FC   WT vs B 14w ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("WT_._14w",colnames(logc))])==0&&sd(logc[i,grep("B_._14w",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("WT_._14w",colnames(logc))]),as.numeric(logc[i,grep("B_._14w",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("WT_._14w",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("B_._14w",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("WT_._14w",colnames(scale_logc)),
                                              grep("B_._14w",colnames(scale_logc)))])), 
                sampleInfo[c(grep("WT_._14w",colnames(scale_logc)),
                             grep("B_._14w",colnames(scale_logc))),]$condition,orthoI = 1) 
vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix,vip)
colnames(mix)
mix$log2FC <- log2(as.numeric(mix$FC))
mix_select=subset(mix,FDR<0.05 & .>1 )
dim(mix_select)
colnames(mix)
res_MF<- mix
#View(res_MF)
res_MF$type <- ifelse(res_MF$FDR< 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')

res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
res_MF$type[which(res_MF$VIP<1)] <- "noSig"
mix_select <- mix_select[order(mix_select$log2FC),]
mix_select$group <- 0
mix_select$group[which(mix_select$log2FC>0)] <- "WildType"
mix_select$group[which(mix_select$log2FC<0)] <- "Blackrice"
Diff_WT_B <- mix_select
Diff_WT_B$Compounds <- rownames(Diff_WT_B)
dim(Diff_WT_B)

bc=Diff_C_B[Diff_C_B$group=="Blackrice",]
cb=Diff_C_B[Diff_C_B$group=="Control",]
wtb=Diff_WT_B[Diff_WT_B$group=="WildType",]
bwt=Diff_WT_B[Diff_WT_B$group=="Blackrice",]
cwt=Diff_C_WT[Diff_C_WT$group=="Control",]
wtc=Diff_C_WT[Diff_C_WT$group=="WildType",]
h22 <- intersect(wtc$Compounds,bc$Compounds)
d22 <- intersect(cwt$Compounds,bwt$Compounds)

h <- intersect(h14,h22)
d <- intersect(d14,d22)
h <- c(h14,h22)
d <- c(d14,d22)
colnames(ALL_sample_data)
health <- ALL_sample_data[,c(2,185,186)][ALL_sample_data[,c(2,185,186)]$Compounds%in%h,]
disease <- ALL_sample_data[,c(2,185,186)][ALL_sample_data[,c(2,185,186)]$Compounds%in%d,]
head(health)
write.table(health,"health.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(disease,"disease.txt", sep = '\t', col.names = NA, quote = FALSE)


phData <- apply(logc[unique(c(health$Compounds,disease$Compounds)),],1,function(x){tapply(x,sampleInfo$group1,mean)}) %>% t() %>% as.data.frame()
colnames(phData)
phData <- phData[,c("WT_14w","C_14w","B_14w","WT_22w","C_22w","B_22w")]
n <- t(scale(t(phData)))
n[n>2]=2
n[n<-2]=-2
pheatmap::pheatmap(n,kmeans_k = 10,fontsize=16,border_color="white")


library(ComplexHeatmap)
library(circlize)

top_anno = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00")),
                       labels =c("WT_14w","C_14w","B_14w","WT_22w","C_22w","B_22w"),
                       labels_gp = gpar(col = "white", fontsize = 14)))

select <- phData[rownames(phData)%in%c("L-Tryptophan","Indole-3-lactic acid","Indole"),]
dim(select)
gene_pos <- which(rownames(phData)  %in% rownames(select))
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,
                                                 labels = rownames(select)))
col_fun = colorRamp2(c(-1, 0, 1), c("cornflowerblue","white","red"))
n=t(scale(t(phData)))
n
p1 <- Heatmap(n,col=col_fun,border = "black",
              show_row_names = FALSE,show_column_names = F,cluster_columns = FALSE,cluster_rows =T,
              column_split = factor(colnames(phData),levels=c("WT_14w","C_14w","B_14w","WT_22w","C_22w","B_22w")),
              #              row_split = Diff$panel,
              top_annotation = top_anno,
              #              left_annotation = left_anno,
              right_annotation = row_anno,row_title = NULL,column_title = NULL,
              heatmap_legend_param = list(
                title = "scale",
                title_position = "leftcenter-rot"))
p1
HM <- draw(p1)  #Show the heatmap
HM
# r.dend <- row_dend(HM)  #Extract row dendrogram   ## 不知道是干嘛的
rcl.list <- row_order(HM) 
write.table(ALL_meta[rownames(n),],"F:/代谢组-蛋白质合作/diff_order_metab.txt", sep = '\t', col.names = NA, quote = FALSE)



# metabolites origin
#metab_origin
Metabolites_Origin_h <- read.csv("MetOrigin_health/02_Origin_Analysis/Metabolites_Origin.csv")
Metabolites_Origin_h$Origin <- gsub("Drug related", "Others",Metabolites_Origin_h$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin_h$Origin <- factor(Metabolites_Origin_h$Origin,levels = c("Host","Microbiota","Co-Metabolism","Food related","Others"  ))
Metabolites_Origin_h$Class <- "Health"
Metabolites_Origin_h$value <- 1
table(Metabolites_Origin_h$Origin)
Metabolites_Origin_d <- read.csv("MetOrigin_disease/02_Origin_Analysis/Metabolites_Origin.csv")
Metabolites_Origin_d$Origin <- gsub("Drug related", "Others",Metabolites_Origin_d$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin_d$Origin <- factor(Metabolites_Origin_d$Origin,levels = c("Host","Microbiota","Co-Metabolism","Food related","Others"  ))
Metabolites_Origin_d$Class <- "Disease"
Metabolites_Origin_d$value <- 1

data <- rbind(Metabolites_Origin_h,Metabolites_Origin_d)
data$Class
data
data <- data%>%dplyr::group_by(Class)%>%dplyr::summarise(n=table(Origin))%>%as.data.frame()
data$Origin <- rep(c("Host", "Microbiota", "Co-Metabolism", "Food related", "Others"),2)
data$Origin <- factor(data$Origin,levels =c("Host", "Microbiota", "Co-Metabolism", "Food related", "Others") )
View(data)
library(dplyr)
library(reshape2)
data$Class <- factor(data$Class,levels = c("Health","Disease"))

ggplot(data, aes(x=Class,y=n,fill=Origin)) +  
  geom_col(position = position_dodge(width = 0.8),width = 0.8) +  
  geom_text(aes(label=n),position = position_dodge(width = 0.8),vjust=-0.)+
  theme_classic()+scale_fill_manual(values = c("#A0B3A1","#4D7587","#85AEC9","#976276","#E1AB8D"))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+xlab("")+ylab("N of Origin")


# Health
Co <- read.csv("MetOrigin_health/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
Micro <- read.csv("MetOrigin_health/03_Function_Analysis/MPEA_Results_Microbiota.csv")
H <- read.csv("MetOrigin_health/03_Function_Analysis/MPEA_Results_Host.csv")
Co$group <- "Co_Metabolism"
Micro$group <- "Micro_Metabolism"
H$group <- "Host_Metabolism"
hp <- rbind(Co,Micro)%>%rbind(.,H)
hp$logPvalue <- -log(hp$Pvalue)
hp <- hp[order(hp$Pvalue,decreasing = T),]
hp$Class <- "Health"

# Disease
Co <- read.csv("MetOrigin_disease/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
Micro <- read.csv("MetOrigin_disease/03_Function_Analysis/MPEA_Results_Microbiota.csv")
H <- read.csv("MetOrigin_disease/03_Function_Analysis/MPEA_Results_Host.csv")
Co$group <- "Co_Metabolism"
Micro$group <- "Micro_Metabolism"
H$group <- "Host_Metabolism"
dp <- rbind(Co,Micro)%>%rbind(.,H)
dp$logPvalue <- -log(dp$Pvalue)
dp <- dp[order(dp$Pvalue,decreasing = T),]
dp$Class <- "Disease"
data <- rbind(hp,dp)
head(data)
data <- subset(data,logPvalue>2.995732)
data$group <- factor(data$group,levels=c("Host_Metabolism","Micro_Metabolism","Co_Metabolism"))
data$Class <- factor(data$Class,levels = c("Health","Disease"))
#View(bp)
ggplot(data[data$Class=="Health",], aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#A0B3A1","#4D7587","#85AEC9"))+theme_classic()+
  ylab("")+xlab("-log Pvalue")+scale_y_discrete(limits=unique(data[data$Class=="Health",]$Name))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 18,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))


ggplot(data[data$Class=="Disease",], aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#A0B3A1","#4D7587","#85AEC9" ))+theme_classic()+
  ylab("")+xlab("-log Pvalue")+scale_y_discrete(limits=rev(unique(data[data$Class=="Disease",]$Name)))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 18,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))



################################# Correlation #####################################
Apc_meta <- read.delim("F:/A_MetaboAnalysis/Black_rice/Black_rice_metab/Apc_meta.txt")
dim(Apc_meta)
colnames(Apc_meta)
B <- Apc_meta[Apc_meta$Class=="B",]
head(B)
intersect(B$Name,wt_all$Compounds)
health <- intersect(Apc_meta$Name,wt_all$Compounds)

a <- setdiff(Apc_meta$Name,name$Compounds)
Apc_meta[Apc_meta$HMDBID%in%a,]
head(phData)
mean_metab <-apply(logc,1,function(x){tapply(x,paste0(sampleInfo$group1,"_",sampleInfo$sex),mean)}) %>% t() %>% as.data.frame()
write.table(mean_metab,"mean_metab.txt", sep = '\t', col.names = NA, quote = FALSE)

diff_metab <- mean_metab[rownames(mean_metab)%in%Apc_meta$Name,]
#diff_metab <- logc[rownames(logc)%in%Apc_meta$Name,]

mean_micro <- apply(otu_order,1,function(x){tapply(x,paste0(Black_riceInfo$group,"_",Black_riceInfo$sex),mean)}) %>% t() %>% as.data.frame()
diff_micro <- mean_micro[rownames(mean_micro)%in%c("Bifidobacterium_pseudolongum",
                                                   "Bacteroides_uniformis",
                                                   "Lactobacillus_johnsonii",
                                                   "Lactobacillus_murinus",
                                                   "Escherichia_coli",
                                                   "Akkermansia_muciniphila",
                                                   "Lactobacillus_reuteri",
                                                   "Turicimonas_muris",
                                                   "Clostridium_innocuum",
                                                   "Klebsiella_pneumonia"),]
diff_micro <- otu_order[rownames(otu_order)%in%c("Bifidobacterium_pseudolongum",
                                                 "Bacteroides_uniformis",
                                                 "Lactobacillus_johnsonii",
                                                 "Lactobacillus_murinus",
                                                 "Escherichia_coli",
                                                 "Akkermansia_muciniphila",
                                                 "Lactobacillus_reuteri",
                                                 "Turicimonas_muris",
                                                 
                                                 "Klebsiella_pneumonia"),]

#dim(diff_metab)
#write.table(diff_metab,"diff_metab.txt", sep = '\t', col.names = NA, quote = FALSE)
#diff_metab <- read.delim("F:/A_MetaboAnalysis/Black_rice/Black_rice_metab/diff_metab.txt", row.names=1)
#intersect(colnames(diff_metab),colnames(diff_micro))
data1 <-  diff_metab[,intersect(colnames(diff_metab),colnames(diff_micro))]
data2 <- diff_micro[,intersect(colnames(diff_metab),colnames(diff_micro))]
#head(diff_metab)
#head(diff_micro)
#dim(otu_order)
#dim(logc)
#diff_metab <-logc[rownames(logc)%in%Apc_meta$Name,]
#diff_micro <- Apc_otu_order[rownames(Apc_otu_order)%in%c("Bifidobacterium_pseudolongum",
#                                                   "Bacteroides_uniformis",
#                                                   "Lactobacillus_johnsonii",
#                                                   "Lactobacillus_murinus",
#                                                   "Escherichia_coli",
#                                                   "Akkermansia_muciniphila",
#                                                   "Lactobacillus_reuteri",
#                                                   "Turicimonas_muris",
#                                                   "Paeniclostridium_sordellii",
#                                                   "Clostridium_innocuum",
#                                                   "Klebsiella_pneumonia"),]


#a <- intersect(colnames(diff_micro),colnames(diff_metab))
#a
#meta <- data.frame(Metabolite_sampleID=a,Microbiome_sampleID=a,
#                   Grouping=paste0(strsplit2(a,"_")[,1],"_",strsplit2(a,"_")[,2]))
#head(meta)
#data1 <- diff_metab[,a]
#head(name)
#data1$Compounds <- rownames(data1)
#data1 <- merge(data1,name,all.x=T,by="Compounds")
#head(data1)
#data2 <- diff_micro[,a]

data1 <- cbind(data1[,grep("C",colnames(data1))],data1[,grep("B",colnames(data1))])
data2 <- cbind(data2[,grep("C",colnames(data2))],data2[,grep("B",colnames(data2))])

#write.table(meta,"meta.txt", sep = '\t', col.names = NA, quote = FALSE)
#write.table(data1,"data_metab.txt", sep = '\t', col.names = NA, quote = FALSE)
#write.table(data2,"data_micro.txt", sep = '\t', col.names = NA, quote = FALSE)
#View(diff_metab)
library(psych)
#diff_metab <- t(scale(t(diff_metab)))
#cor <-corr.test(t(diff_metab[,a]),t(diff_micro[,a]), method = "spearman",adjust= "none")
cor <-corr.test(t(data1),t(data2), method = "spearman",adjust= "none")
#colnames(diff_metab[,a])
#colnames(diff_micro[,a])
#cor <-corr.test(diff_microb, diff_metab, method = "pearson",adjust= "none")
## 提取相关性、p值

cmt <-cor$r
head(cmt)
c=melt(cmt)
pmt <- cor$p
#p=melt(pmt)
#cp <- cbind(c,p[,3])

## 输出相关系数表格,第一行为代谢物信息，第一列为物种信息
cmt.out<-cbind(rownames(cmt),cmt)
#write.table(cmt.out,file= "cor.txt",sep= "t",row.names=F)
## 输出p值表格，第一行为代谢物信息，第一列为物种信息
pmt.out<-cbind(rownames(pmt),pmt)
#write.table(pmt.out,file= "pvalue.txt",sep= "t",row.names=F)
## 第一列为物种名，第二列为代谢物名，第三、第四列对应显示相关系数与p值
df <-melt(cmt,value.name= "cor")
df$pvalue <- as.vector(pmt)
head(df)
df$fdr <- p.adjust(df$pvalue,method = "fdr")
df <- subset(df,abs(cor)>0.6& fdr<0.05)
#View(df)
#write.table(df,file= "cor-p.txt",sep= "t")

if(!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <- '**'
  smt <- pmt > 0.01& pmt < 0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else{
  pmt <- F
}
pmt
mycol <-  colorRampPalette(c("#336699", "white", "#CC3333"))(200)
#data <- t(cmt)
#data[,109]
#apply(data,1,which.max)
#a=unique(Apc_meta[Apc_meta$Class%in%c("B","W","P","A","F"),]$Name)
#length(a)
#dim(cmt)
#a=intersect(a,rownames(cmt))
#data <- cmt[a,]
library(ComplexHeatmap)
library(circlize)
select <- data[rownames(data)%in%c(health,"L-Histidine","L-Tryptophan","Indole-3-lactic acid",
                                   "Estriol","Biopterin","D-Tagatose","4-Hydroxybenzoic Acid",
                                   "14,15-DHET","Riboflavin","Phosphocholine","Malonicacid"),]


select <- data[rownames(data)%in%c("L-Tryptophan","Indole-3-lactic acid",
                                   "Indole"),]
select
#"L-Alanine","L-Glutamic Acid","D-Galactose","Myoinositol",
#"ADP-ribose","Adenosine","Hypoxanthine","Xanthine"
select <- select[!duplicated(rownames(select)),]
dim(select)
gene_pos <- which(rownames(data)  %in% rownames(select))
gene_pos <- which(unique(rownames(data))%in% unique(rownames(select)))
col_fun = colorRamp2(c(-0.6,0,0.6), c("cornflowerblue","white","red"))
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,
                                                 labels = unique(rownames(select))))
colnames(data)
p1 <- Heatmap(cmt[c("L-Tryptophan","Indole-3-lactic acid","Indole"),],col=col_fun,border = "black",
              show_row_names = T,show_column_names =T,
              cluster_columns = F,cluster_rows = T,
              #left_annotation = left_anno,
              right_annotation = row_anno,
              row_title = NULL,column_title = NULL,
              heatmap_legend_param = list(title = "scale",
                                          title_position = "leftcenter-rot"))
p1


pheatmap::pheatmap(cmt[c("L-Tryptophan","Indole-3-lactic acid","Indole"),],fontsize = 20,legend_breaks = c(-1,-0.5,0,0.5,1),display_numbers = t(pmt),
                   color=mycol)

diff_micro[c("Escherichia_coli","Bacteroides_uniformis"),]
diff_metab[c("L-Tryptophan","Indole","Indole-3-lactic acid"),]



ggdata <- t(rbind(data1,data2))%>%as.data.frame(.)
#dim(data)
#dim(logh[,-78])
colnames(ggdata)
View(ggdata)
ggplot(data=ggdata[44:73,], aes(x=`4-hydroxyphenylglyoxylate`, y=`Bacteroides_uniformis`))+geom_point(color="black",size=3)+
  stat_smooth(method="lm",se=T)+
  stat_cor(data=ggdata[44:73,], method = "spearman")+theme_classic()+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))
rownames(ggdata)
ggplot(data=ggdata[30:73,], aes(x=`4-hydroxyphenylglyoxylate`, y=`Akkermansia_muciniphila`))+geom_point(color="black",size=3)+
  stat_smooth(method="lm",se=T)+
  stat_cor(data=ggdata[30:73,], method = "spearman")+theme_classic()+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))
View(ALL_sample_data)


remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}

data <- logc[rownames(logc)%in%c("L-Tryptophan","Indole",
                                 "Indole-3-lactic acid"),]
data$metab <- rownames(data)
ggdata <- melt(data)

ggdata
meta <- strsplit2(ggdata$variable,"_")%>%as.data.frame()
ggdata$group <- paste0(meta[,1],"_",meta[,3])%>%
  factor(.,levels = c("WT_14w","C_14w","B_14w",
                      "WT_22w","C_22w","B_22w"))
data <- ggdata[ggdata$metab=="L-Tryptophan",]
data1 <- data %>% 
  dplyr::group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()   
ggplot(na.omit(data1),aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  ylab("Relative concentration")+ggtitle("L-Tryptophan")+ 
  theme(plot.title = element_text(size = 24, face = "bold"))

data <- ggdata[ggdata$metab=="Indole",]
data1 <- data %>% 
  dplyr::group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()   
ggplot(na.omit(data1),aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  ylab("Relative concentration")+ggtitle("Indole")+ 
  theme(plot.title = element_text(size = 24, face = "bold"))

data <- ggdata[ggdata$metab=="Indole-3-lactic acid",]
data1 <- data %>% 
  dplyr::group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()   
ggplot(na.omit(data1),aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  ylab("Relative concentration")+ggtitle("Indole-3-lactic acid")+ 
  theme(plot.title = element_text(size = 24, face = "bold"))



######################### 中介效应 分析 ##############################
data1 <-  diff_metab[,intersect(colnames(diff_metab),colnames(diff_micro))]
data2 <- diff_micro[,intersect(colnames(diff_metab),colnames(diff_micro))]
data <- rbind(data1,data2)
df <- data[rownames(data)%in%c("L-Tryptophan",
                               "Indole",
                               "Indole-3-lactic acid","Bacteroides_uniformis"),]%>%t()%>%as.data.frame()

df$group <- 0
df$num <- 1:length(rownames(df))
View(df)
df$group[1:14] <- 1
df$group[74:88] <- 1
df <- df[grep("WT",rownames(df),invert = T),]
df

df$group<- as.factor(df$group)
#建立线性回归，我的X是CX3CL1, M是Meta52, Y是Group（二分类变量）
a<- lm(Bacteroides_uniformis ~ group,df) #lm(M~X,df)
b<- lm(`L-Tryptophan`~group+Bacteroides_uniformis, df) #glm(Y~X+M)
#install.packages("mediation")
library(mediation)
set.seed(123) #保证结果可以复现
result = mediate(a,b,treat="group",mediator = "Bacteroides_uniformis",boot = T)#默认1000次抽样
#这里需要注意变量的名称，如果名称比较特殊，有中文或者符号很容易出错，建议提前处理变量名
summary(result)

data <- rbind(data1,data2)%>%t()%>%as.data.frame()
data <- data[,colnames(data)%in%c("L-Tryptophan",
                                  "Indole",
                                  "Indole-3-lactic acid","Bacteroides_uniformis")]


ggplot(data=data, aes(x=`Bacteroides_uniformis`, y=`Indole`))+
  geom_point(color="black",size=5)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=data, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))


ggplot(data=data, aes(x=`Bacteroides_uniformis`, y=`Indole-3-lactic acid`))+
  geom_point(color="black",size=5)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=data, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))


data$num <- 1:length(rownames(data))
View(data)
data <- data.frame(Bacteroides_uniformis=rep(data$Bacteroides_uniformis,3),Tryptophan=c(data$Indole,data$`Indole-3-lactic acid`,data$`L-Tryptophan`))
data$group <- c(rep("Indole",20),rep("Indole-3-lactic acid",20),rep("L-Tryptophan",20))
ggscatter(data, x = "Bacteroides_uniformis", y = "Tryptophan",
          add = "reg.line",
          conf.int = TRUE,                         
          color = "group", palette = "jco"         
) +stat_cor(method = "spearman")   
