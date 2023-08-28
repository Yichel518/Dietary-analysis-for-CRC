setwd("F:/A_MetaboAnalysis/Black_rice/Black_rice_metab/AOM_results/")
#################################### Apc BR #########################################
library(readxl)
library(tidyverse)
library(limma)
ALL_sample_data <- read_excel("F:/A_MetaboAnalysis/Black_rice/Black_rice_metab/AOM_results/1.Data_Assess/all_group/ALL_sample_data.xlsx")
colnames(ALL_sample_data)
ALL_meta <- data.frame(Compounds=ALL_sample_data$Compounds,ALL_sample_data[,12:186])


colnames(ALL_meta)
ALL_meta[1:5,1:5]
rownames(ALL_meta)<- ALL_meta$Compounds
ALL_meta$Compounds <- NULL
head(ALL_meta)
name <- ALL_sample_data[,c(2,209)]
ALL_meta <- ALL_meta[,order(colnames(ALL_meta))]
sampleInfo <- data.frame(strsplit2(colnames(ALL_meta),"_"))
colnames(sampleInfo) <- c("condition","sex","stage","replicate")
#sampleInfo$group1 <- paste0(sampleInfo$condition,"_",sampleInfo$sex)
sampleInfo$group <- paste0(sampleInfo$condition,"_",sampleInfo$stage)
#sampleInfo$group2 <- paste0(sampleInfo$condition,"_",sampleInfo$sex,"_",sampleInfo$stage)
rownames(sampleInfo) <- colnames(ALL_meta)
sampleInfo$group<- factor(sampleInfo$group,levels = c("WT_1","WT_3","C_1","C_3","B_1","B_3",
                                                      "W_1","W_3","P_1","P_3"))

dim(sampleInfo)
colnames(ALL_meta)
logc <- log(ALL_meta+1)
head(logc)
rownames(ALL_meta)[1:5]
dim(logc)
scale_logc <- t(scale(t(logc)))
library(PCAtools)
dim(logc)
dim(sampleInfo)
pca <- PCAtools::pca(as.data.frame(logc), metadata =sampleInfo)
biplot(pca, x = 'PC1', y = 'PC2') #可以看到这两个成分对样品的解释度
#screen(p) #可以看到所有成分对样品的解释
##将pca与关联样本信息表 
pca$variance
pca_rotated_plus <- cbind(pca$rotated,sampleInfo)
library(ggsci)
pca_rotated_plus$group <- factor(pca_rotated_plus$group,levels = c("WT_1","WT_3","C_1","C_3","B_1","B_3",
ggplot(pca_rotated_plus,aes(x = PC1 , y = PC2))+
  geom_point(size = 8,aes(fill = group,color = group)) +
  stat_ellipse(aes(color = group,fill=group),linetype = 'dashed',size = 1, show.legend = FALSE) + #添加分组椭圆
  labs(x = 'PC1 (22.48%)',y = 'PC2 (13.31%)') + 
  #    scale_shape_manual(values = c(21,21,21,21,21,22,22,22,22,22))+ 
  scale_fill_manual(values = c("#83BFE3","#125AAF","#FFDD93","#F5C31E","#0ECBC5",
                               "#0293A1","#E56A54","#F73E00","#D5EFA7","#56D75E"))+
  scale_color_manual(values = c("#83BFE3","#125AAF","#FFDD93","#F5C31E","#0ECBC5",
                                "#0293A1","#E56A54","#F73E00","#D5EFA7","#56D75E"))+
  theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))+
  theme(legend.title=element_blank())

########################### pvalue FC C vs WT Batch 1  ################################
Pvalue<-c(rep(0,nrow(logc))) 
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._1_.",colnames(logc))])==0&&sd(logc[i,grep("WT_._1_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._1_.",colnames(logc))]),as.numeric(logc[i,grep("WT_._1_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._1_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("WT_._1_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)
head(mix)
mix$log2FC <- log2(as.numeric(mix$FC))
#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("C_._1_.",colnames(scale_logc)),
                                              grep("WT_._1_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("C_._1_.",colnames(scale_logc)),
                             grep("WT_._1_.",colnames(scale_logc))),]$condition,orthoI = 1) 
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

######################## pvalue FC   C vs B ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._1_.",colnames(logc))])==0&&sd(logc[i,grep("B_._1_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._1_.",colnames(logc))]),as.numeric(logc[i,grep("B_._1_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._1_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("B_._1_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("C_._1_.",colnames(scale_logc)),
                                              grep("B_._1_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("C_._1_.",colnames(scale_logc)),
                             grep("B_._1_.",colnames(scale_logc))),]$condition,orthoI = 1) 
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
######################## pvalue FC   C vs P ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._1_.",colnames(logc))])==0&&sd(logc[i,grep("P_._1_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._1_.",colnames(logc))]),as.numeric(logc[i,grep("P_._1_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._1_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("P_._1_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("C_._1_.",colnames(scale_logc)),
                                              grep("P_._1_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("C_._1_.",colnames(scale_logc)),
                             grep("P_._1_.",colnames(scale_logc))),]$condition,orthoI = 1) 
vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
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
mix_select$group[which(mix_select$log2FC<0)] <- "Polishrice"
Diff_C_P <- mix_select
Diff_C_P$Compounds <- rownames(Diff_C_P)
dim(Diff_C_P)
######################## pvalue FC   C vs W ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._1_.",colnames(logc))])==0&&sd(logc[i,grep("W_._1_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._1_.",colnames(logc))]),as.numeric(logc[i,grep("W_._1_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._1_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("W_._1_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("C_._1_.",colnames(scale_logc)),
                                              grep("W_._1_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("C_._1_.",colnames(scale_logc)),
                             grep("W_._1_.",colnames(scale_logc))),]$condition,orthoI = 1) 

vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
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
mix_select$group[which(mix_select$log2FC<0)] <- "Brownrice"
Diff_C_W <- mix_select
Diff_C_W$Compounds <- rownames(Diff_C_W)
dim(Diff_C_W)
head(Diff_C_B)
head(Diff_C_P)
head(Diff_C_W)

a_1=Diff_C_B[Diff_C_B$group=="Blackrice",]
dim(a)
#View(a)
b_1=Diff_C_P[Diff_C_P$group=="Polishrice",]
dim(b_1)
c_1=Diff_C_W[Diff_C_W$group=="Brownrice",]
dim(c_1)
d_1=Diff_C_WT[Diff_C_WT$group=="WildType",]
dim(d_1)
ggdata = list(rownames(a_1),rownames(b_1),rownames(c_1),rownames(d_1))

names(ggdata) <- c("B","P","W","WT")

#BiocManager::install("ggvenn")
library(ggvenn)
p=ggvenn(ggdata,       
         show_percentage = F,show_elements=F,text_size = 8,
         stroke_color = "Black",
         fill_color = c("#85AEC9","#A0B3A1","#F8C42C","#E56A54"),
         set_name_color = c("Black","Black","Black","Black"),set_name_size = 8) 
p
dim(name)
B_1=setdiff(rownames(a_1),c(rownames(b_1),rownames(c_1)))
length(B_1)
P_1=setdiff(rownames(b_1),c(rownames(a_1),rownames(c_1)))
length(P_1)
W_1=setdiff(rownames(c_1),c(rownames(a_1),rownames(b_1)))
length(W_1)
S_1=intersect(rownames(a_1),rownames(b_1))%>%intersect(.,rownames(c_1))
length(S_1)
B_specific_1 <- name[name$Compounds%in%B_1,]
P_specific_1 <- name[name$Compounds%in%P_1,]
W_specific_1 <- name[name$Compounds%in%W_1,]
S_specific_1 <- name[name$Compounds%in%S_1,]

write.table(B_specific_1,'B_specific_1.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(P_specific_1,'P_specific_1.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(W_specific_1,'W_specific_1.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(S_specific_1,'S_specific_1.txt', sep = '\t', col.names = NA, quote = FALSE)

c("Indole-3-lactic acid","Cabergoline")
table(class_br$`Class I`)
class_br <- ALL_sample_data[ALL_sample_data$Compounds%in%a_1$Compounds,][,c(2,4,6)]%>%as.data.frame()
View(class_br)
class_br <- as.data.frame(table(class_br$`Class II`))
class_br <- class_br[order(class_br$Freq),]

ggplot(data=class_br,aes(x=Freq,y=Var1)) +
  geom_bar(data=class_br,aes(x=Freq,y=Var1), stat = "identity")+
  ylab("Metablites Class")+xlab("Hits Count")+scale_y_discrete(limits=class_br$Var1)+ theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))


class_pr <- ALL_sample_data[ALL_sample_data$Compounds%in%P,][,c(2,4,6)]%>%as.data.frame()
class_pr <- as.data.frame(table(class_pr$`Class II`))
class_pr <- class_pr[order(class_pr$Freq),]
ggplot(data=class_pr,aes(x=Freq,y=Var1)) +
  geom_bar(data=class_pr,aes(x=Freq,y=Var1), stat = "identity")+
  ylab("Metablites Class")+xlab("Hits Count")+scale_y_discrete(limits=class_pr$Var1)+ theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))+scale_x_continuous(limits=c(0,50), breaks=seq(0,50,10))

class_wr <- ALL_sample_data[ALL_sample_data$Compounds%in%W,][,c(2,4,6)]%>%as.data.frame()
class_wr <- as.data.frame(table(class_wr$`Class II`))
class_wr <- class_wr[order(class_wr$Freq),]

ggplot(data=class_wr,aes(x=Freq,y=Var1)) +
  geom_bar(data=class_wr,aes(x=Freq,y=Var1), stat = "identity")+
  ylab("Metablites Class")+xlab("Hits Count")+scale_y_discrete(limits=class_wr$Var1)+ theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))+scale_x_continuous(limits=c(0,50), breaks=seq(0,50,10))











######################## pvalue FC   B vs W ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("B_._1_.",colnames(logc))])==0&&sd(logc[i,grep("W_._1_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("B_._1_.",colnames(logc))]),as.numeric(logc[i,grep("W_._1_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("B_._1_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("W_._1_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("B_._1_.",colnames(scale_logc)),
                                              grep("W_._1_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("B_._1_.",colnames(scale_logc)),
                             grep("W_._1_.",colnames(scale_logc))),]$condition,orthoI = 1) 

vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
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
mix_select$group[which(mix_select$log2FC>0)] <- "Blackrice"
mix_select$group[which(mix_select$log2FC<0)] <- "Brownrice"
Diff_B_W <- mix_select
Diff_B_W$Compounds <- rownames(Diff_B_W)
dim(Diff_B_W)

######################## pvalue FC   P vs W ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("P_._1_.",colnames(logc))])==0&&sd(logc[i,grep("W_._1_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("P_._1_.",colnames(logc))]),as.numeric(logc[i,grep("W_._1_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("P_._1_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("W_._1_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("P_._1_.",colnames(scale_logc)),
                                              grep("W_._1_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("P_._1_.",colnames(scale_logc)),
                             grep("W_._1_.",colnames(scale_logc))),]$condition,orthoI = 1) 

vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
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
mix_select$group[which(mix_select$log2FC>0)] <- "polisehrice"
mix_select$group[which(mix_select$log2FC<0)] <- "Brownrice"
Diff_P_W <- mix_select
Diff_P_W$Compounds <- rownames(Diff_P_W)
dim(Diff_P_W)

######################## pvalue FC   B vs P ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("B_._1_.",colnames(logc))])==0&&sd(logc[i,grep("P_._1_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("B_._1_.",colnames(logc))]),as.numeric(logc[i,grep("P_._1_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("B_._1_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("P_._1_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("B_._1_.",colnames(scale_logc)),
                                              grep("P_._1_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("B_._1_.",colnames(scale_logc)),
                             grep("P_._1_.",colnames(scale_logc))),]$condition,orthoI = 1) 

vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
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
mix_select$group[which(mix_select$log2FC>0)] <- "Blackrice"
mix_select$group[which(mix_select$log2FC<0)] <- "Polishrice"
Diff_B_P <- mix_select
Diff_B_P$Compounds <- rownames(Diff_B_P)
dim(Diff_B_P)
#View(Diff_C_WT)
#View(Diff_C_B)
#View(Diff_C_P)
#View(Diff_C_W)
#View(Diff_B_W)
#View(Diff_P_W)
Diff_C_WT$metab <- rownames(Diff_C_WT)
Diff_C_WT$panel <- "CvsWT"
Diff_C_B$metab <- rownames(Diff_C_B)
Diff_C_B$panel <- "CvsB"
Diff_C_P$metab <- rownames(Diff_C_P)
Diff_C_P$panel <- "CvsP"
Diff_C_W$metab <- rownames(Diff_C_W)
Diff_C_W$panel <- "CvsW"
Diff_B_W$metab <- rownames(Diff_B_W)
Diff_B_W$panel <- "BvsW"
Diff_P_W$metab <- rownames(Diff_P_W)
Diff_P_W$panel <- "PvsW"
write.table(Diff_C_WT,"Diff_C_WT_1.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(Diff_C_B,"Diff_C_B_1.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(Diff_C_P,"Diff_C_P_1.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(Diff_C_W,"Diff_C_W_1.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(Diff_B_W,"Diff_B_W_1.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(Diff_P_W,"Diff_P_W_1.txt", sep = '\t', col.names = NA, quote = FALSE)
Diff_1 <- rbind(Diff_C_WT,Diff_C_B)%>%rbind(.,Diff_C_P)%>%rbind(.,Diff_C_W)%>%rbind(.,Diff_B_W)%>%rbind(.,Diff_P_W)
Diff_1$panel <- factor(Diff_1$panel,levels = c("CvsWT","CvsB","CvsW","CvsP","PvsW","BvsW"))
write.table(Diff_1,"Diff_1.txt", sep = '\t', col.names = NA, quote = FALSE)



########################### pvalue FC C vs WT Batch 3  ################################
Pvalue<-c(rep(0,nrow(logc))) 
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._3_.",colnames(logc))])==0&&sd(logc[i,grep("WT_._3_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._3_.",colnames(logc))]),as.numeric(logc[i,grep("WT_._3_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._3_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("WT_._3_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)
head(mix)
mix$log2FC <- log2(as.numeric(mix$FC))
#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("C_._3_.",colnames(scale_logc)),
                                              grep("WT_._3_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("C_._3_.",colnames(scale_logc)),
                             grep("WT_._3_.",colnames(scale_logc))),]$condition,orthoI = 1) 
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

######################## pvalue FC   C vs B ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._3_.",colnames(logc))])==0&&sd(logc[i,grep("B_._3_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._3_.",colnames(logc))]),as.numeric(logc[i,grep("B_._3_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._3_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("B_._3_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("C_._3_.",colnames(scale_logc)),
                                              grep("B_._3_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("C_._3_.",colnames(scale_logc)),
                             grep("B_._3_.",colnames(scale_logc))),]$condition,orthoI = 1) 
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
######################## pvalue FC   C vs P ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._3_.",colnames(logc))])==0&&sd(logc[i,grep("P_._3_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._3_.",colnames(logc))]),as.numeric(logc[i,grep("P_._3_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._3_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("P_._3_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("C_._3_.",colnames(scale_logc)),
                                              grep("P_._3_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("C_._3_.",colnames(scale_logc)),
                             grep("P_._3_.",colnames(scale_logc))),]$condition,orthoI = 1) 
vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
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
mix_select$group[which(mix_select$log2FC<0)] <- "Polishrice"
Diff_C_P <- mix_select
Diff_C_P$Compounds <- rownames(Diff_C_P)
dim(Diff_C_P)
######################## pvalue FC   C vs W ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("C_._3_.",colnames(logc))])==0&&sd(logc[i,grep("W_._3_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("C_._3_.",colnames(logc))]),as.numeric(logc[i,grep("W_._3_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("C_._3_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("W_._3_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("C_._3_.",colnames(scale_logc)),
                                              grep("W_._3_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("C_._3_.",colnames(scale_logc)),
                             grep("W_._3_.",colnames(scale_logc))),]$condition,orthoI = 1) 

vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
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
mix_select$group[which(mix_select$log2FC<0)] <- "Brownrice"
Diff_C_W <- mix_select
Diff_C_W$Compounds <- rownames(Diff_C_W)
dim(Diff_C_W)
head(Diff_C_B)
head(Diff_C_P)
head(Diff_C_W)

a_3=Diff_C_B[Diff_C_B$group=="Blackrice",]
dim(a_3)
#View(a)
b_3=Diff_C_P[Diff_C_P$group=="Polishrice",]
c_3=Diff_C_W[Diff_C_W$group=="Brownrice",]
d_3=Diff_C_WT[Diff_C_WT$group=="WildType",]
ggdata = list(rownames(a_3),rownames(b_3),rownames(c_3),rownames(d_3))
names(ggdata) <- c("B","P","W","WT")
ggdata
#BiocManager::install("ggvenn")
library(ggvenn)
p=ggvenn(ggdata,       
         show_percentage = F,show_elements=F,text_size = 8,
         stroke_color = "Black",
         fill_color = c("#85AEC9","#A0B3A1","#F8C42C","#E56A54"),
         set_name_color = c("Black","Black","Black","Black"),set_name_size = 8) 
p

B_3=setdiff(rownames(a_3),c(rownames(b_3),rownames(c_3)))
length(B_3)
P_3=setdiff(rownames(b_3),c(rownames(a_3),rownames(c_3)))
length(P_3)
W_3=setdiff(rownames(c_3),c(rownames(a_3),rownames(b_3)))
length(W_3)
S_3=intersect(rownames(a_3),rownames(b_3))%>%intersect(.,rownames(c_3))
length(S_3)
B_specific_3 <- name[name$Compounds%in%B_3,]
P_specific_3 <- name[name$Compounds%in%P_3,]
W_specific_3 <- name[name$Compounds%in%W_3,]
S_specific_3 <- name[name$Compounds%in%S_3,]
write.table(B_specific_3,'B_specific_3.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(P_specific_3,'P_specific_3.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(W_specific_3,'W_specific_3.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(S_specific_3,'S_specific_3.txt', sep = '\t', col.names = NA, quote = FALSE)

table(class_br$`Class I`)
class_br <- ALL_sample_data[ALL_sample_data$Compounds%in%a$Compounds,][,c(2,4,6)]%>%as.data.frame()
View(class_br)
class_br <- as.data.frame(table(class_br$`Class II`))
class_br <- class_br[order(class_br$Freq),]

ggplot(data=class_br,aes(x=Freq,y=Var1)) +
  geom_bar(data=class_br,aes(x=Freq,y=Var1), stat = "identity")+
  ylab("Metablites Class")+xlab("Hits Count")+scale_y_discrete(limits=class_br$Var1)+ theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))


class_pr <- ALL_sample_data[ALL_sample_data$Compounds%in%P,][,c(2,4,6)]%>%as.data.frame()
class_pr <- as.data.frame(table(class_pr$`Class II`))
class_pr <- class_pr[order(class_pr$Freq),]
ggplot(data=class_pr,aes(x=Freq,y=Var1)) +
  geom_bar(data=class_pr,aes(x=Freq,y=Var1), stat = "identity")+
  ylab("Metablites Class")+xlab("Hits Count")+scale_y_discrete(limits=class_pr$Var1)+ theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))+scale_x_continuous(limits=c(0,50), breaks=seq(0,50,10))

class_wr <- ALL_sample_data[ALL_sample_data$Compounds%in%W,][,c(2,4,6)]%>%as.data.frame()
class_wr <- as.data.frame(table(class_wr$`Class II`))
class_wr <- class_wr[order(class_wr$Freq),]

ggplot(data=class_wr,aes(x=Freq,y=Var1)) +
  geom_bar(data=class_wr,aes(x=Freq,y=Var1), stat = "identity")+
  ylab("Metablites Class")+xlab("Hits Count")+scale_y_discrete(limits=class_wr$Var1)+ theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))+scale_x_continuous(limits=c(0,50), breaks=seq(0,50,10))











######################## pvalue FC   B vs W ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("B_._3_.",colnames(logc))])==0&&sd(logc[i,grep("W_._3_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("B_._3_.",colnames(logc))]),as.numeric(logc[i,grep("W_._3_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("B_._3_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("W_._3_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("B_._3_.",colnames(scale_logc)),
                                              grep("W_._3_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("B_._3_.",colnames(scale_logc)),
                             grep("W_._3_.",colnames(scale_logc))),]$condition,orthoI = 1) 

vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
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
mix_select$group[which(mix_select$log2FC>0)] <- "Blackrice"
mix_select$group[which(mix_select$log2FC<0)] <- "Brownrice"
Diff_B_W <- mix_select
Diff_B_W$Compounds <- rownames(Diff_B_W)
dim(Diff_B_W)

######################## pvalue FC   P vs W ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("P_._3_.",colnames(logc))])==0&&sd(logc[i,grep("W_._3_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("P_._3_.",colnames(logc))]),as.numeric(logc[i,grep("W_._3_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("P_._3_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("W_._3_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("P_._3_.",colnames(scale_logc)),
                                              grep("W_._3_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("P_._3_.",colnames(scale_logc)),
                             grep("W_._3_.",colnames(scale_logc))),]$condition,orthoI = 1) 

vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
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
mix_select$group[which(mix_select$log2FC>0)] <- "polisehrice"
mix_select$group[which(mix_select$log2FC<0)] <- "Brownrice"
Diff_P_W <- mix_select
Diff_P_W$Compounds <- rownames(Diff_P_W)
dim(Diff_P_W)

######################## pvalue FC   B vs P ################################
Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,grep("B_._3_.",colnames(logc))])==0&&sd(logc[i,grep("P_._3_.",colnames(logc))])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,grep("B_._3_.",colnames(logc))]),as.numeric(logc[i,grep("P_._3_.",colnames(logc))]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,grep("B_._3_.",colnames(ALL_meta))]))/mean(as.numeric(ALL_meta[i,grep("P_._3_.",colnames(ALL_meta))])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

#View(mix)
# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc[,c(grep("B_._3_.",colnames(scale_logc)),
                                              grep("P_._3_.",colnames(scale_logc)))])), 
                sampleInfo[c(grep("B_._3_.",colnames(scale_logc)),
                             grep("P_._3_.",colnames(scale_logc))),]$condition,orthoI = 1) 

vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)
head(mix)
mix <- cbind(mix[rownames(vip),],vip)
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
mix_select$group[which(mix_select$log2FC>0)] <- "Blackrice"
mix_select$group[which(mix_select$log2FC<0)] <- "Polishrice"
Diff_B_P <- mix_select
Diff_B_P$Compounds <- rownames(Diff_B_P)
dim(Diff_B_P)
#View(Diff_C_WT)
#View(Diff_C_B)
#View(Diff_C_P)
#View(Diff_C_W)
#View(Diff_B_W)
#View(Diff_P_W)
Diff_C_WT$metab <- rownames(Diff_C_WT)
Diff_C_WT$panel <- "CvsWT"
Diff_C_B$metab <- rownames(Diff_C_B)
Diff_C_B$panel <- "CvsB"
Diff_C_P$metab <- rownames(Diff_C_P)
Diff_C_P$panel <- "CvsP"
Diff_C_W$metab <- rownames(Diff_C_W)
Diff_C_W$panel <- "CvsW"
Diff_B_W$metab <- rownames(Diff_B_W)
Diff_B_W$panel <- "BvsW"
Diff_P_W$metab <- rownames(Diff_P_W)
Diff_P_W$panel <- "PvsW"
write.table(Diff_C_WT,"Diff_C_WT_3.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(Diff_C_B,"Diff_C_B_3.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(Diff_C_P,"Diff_C_P_3.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(Diff_C_W,"Diff_C_W_3.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(Diff_B_W,"Diff_B_W_3.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(Diff_P_W,"Diff_P_W_3.txt", sep = '\t', col.names = NA, quote = FALSE)
Diff_3 <- rbind(Diff_C_WT,Diff_C_B)%>%rbind(.,Diff_C_P)%>%rbind(.,Diff_C_W)%>%rbind(.,Diff_B_W)%>%rbind(.,Diff_P_W)
Diff_3$panel <- factor(Diff_3$panel,levels = c("CvsWT","CvsB","CvsW","CvsP","PvsW","BvsW"))
write.table(Diff_3,"Diff_3.txt", sep = '\t', col.names = NA, quote = FALSE)


































#BR
br_inte <- name[name$Compounds%in%intersect(B_specific_1$Compounds,
                                            B_specific_3$Compounds),]
br_all=name[name$Compounds%in%c(B_specific_1$Compounds,B_specific_3$Compounds),]
br_1 <- name[name$Compounds%in%setdiff(B_specific_1$Compounds,B_specific_3$Compounds),]
br_3 <- name[name$Compounds%in%setdiff(B_specific_3$Compounds,B_specific_1$Compounds),]
write.table(br_all,"br_all.txt", sep = '\t', col.names = NA, quote = FALSE)
length(br_all$Compounds)
length(br_inte$Compounds)
length(br_1$Compounds)
length(br_3$Compounds)

# WR
wr_inte <- name[name$Compounds%in%intersect(W_specific_1$Compounds,
                                            W_specific_3$Compounds),]
wr_all=name[name$Compounds%in%c(W_specific_1$Compounds,W_specific_3$Compounds),]
wr_1 <- name[name$Compounds%in%setdiff(W_specific_1$Compounds,W_specific_3$Compounds),]
wr_3 <- name[name$Compounds%in%setdiff(W_specific_3$Compounds,W_specific_1$Compounds),]
write.table(wr_all,"wr_all.txt", sep = '\t', col.names = NA, quote = FALSE)
length(wr_all$Compounds)
length(wr_inte$Compounds)
length(wr_1$Compounds)
length(wr_3$Compounds)

# PR
pr_inte <- name[name$Compounds%in%intersect(P_specific_1$Compounds,
                                            P_specific_3$Compounds),]
pr_all=name[name$Compounds%in%c(P_specific_1$Compounds,P_specific_3$Compounds),]
pr_1 <- name[name$Compounds%in%setdiff(P_specific_1$Compounds,P_specific_3$Compounds),]
pr_3 <- name[name$Compounds%in%setdiff(P_specific_3$Compounds,P_specific_1$Compounds),]
write.table(pr_all,"pr_all.txt", sep = '\t', col.names = NA, quote = FALSE)
length(pr_all$Compounds)
length(pr_inte$Compounds)
length(pr_1$Compounds)
length(pr_3$Compounds)

# shared 
## 1
b1 <- Diff_1[Diff_1$panel%in%c("CvsB"),]%>%.[.$group=="Blackrice",]
b1 <- name[name$Compounds%in%b1$Compounds,]
w1 <- Diff_1[Diff_1$panel%in%c("CvsW"),]%>%.[.$group=="Brownrice",]
w1 <- name[name$Compounds%in%w1$Compounds,]
p1 <- Diff_1[Diff_1$panel%in%c("CvsP"),]%>%.[.$group=="Polishrice",]
p1 <- name[name$Compounds%in%p1$Compounds,]
s1 <- intersect(b1$Compounds,w1$Compounds)%>%intersect(.,p1$Compounds)
s1 <- name[name$Compounds%in%s1,]
dim(s1)

# shared 
## 2
b3 <- Diff_3[Diff_3$panel%in%c("CvsB"),]%>%.[.$group=="Blackrice",]
b3 <- name[name$Compounds%in%b3$Compounds,]
w3 <- Diff_3[Diff_3$panel%in%c("CvsW"),]%>%.[.$group=="Brownrice",]
w3 <- name[name$Compounds%in%w3$Compounds,]
p3 <- Diff_3[Diff_3$panel%in%c("CvsP"),]%>%.[.$group=="Polishrice",]
p3 <- name[name$Compounds%in%p3$Compounds,]
s3 <- intersect(b3$Compounds,w3$Compounds)%>%intersect(.,p3$Compounds)
s3 <- name[name$Compounds%in%s3,]
dim(s3)
s_inte <- name[name$Compounds%in%intersect(s1$Compounds,s3$Compounds),]
s_all <- name[name$Compounds%in%c(s1$Compounds,s3$Compounds),]
s_1 <- name[name$Compounds%in%setdiff(s1$Compounds,s3$Compounds),]
s_3 <- name[name$Compounds%in%setdiff(s3$Compounds,s1$Compounds),]
write.table(s_all,"s_all.txt", sep = '\t', col.names = NA, quote = FALSE)
length(s_all$Compounds)
length(s_inte$Compounds)
length(s_1$Compounds)
length(s_3$Compounds)

# 花青素（a）
View(Diff_1[Diff_1$panel%in%c("BvsW"),])
a1 <- Diff_1[Diff_1$panel%in%c("BvsW"),]%>%.[.$group=="Blackrice",]
a3 <- Diff_3[Diff_3$panel%in%c("BvsW"),]%>%.[.$group=="Blackrice",]
length(a1$Compounds)
length(a3$Compounds)
a_inte <- name[name$Compounds%in%intersect(a1$Compounds,a3$Compounds),]
a_all <- name[name$Compounds%in%c(a1$Compounds,a3$Compounds),]
a_1 <- name[name$Compounds%in%setdiff(a1$Compounds,a3$Compounds),]
a_3 <- name[name$Compounds%in%setdiff(a3$Compounds,a1$Compounds),]

write.table(a_all,"a_all.txt", sep = '\t', col.names = NA, quote = FALSE)
length(a_all$Compounds)
length(a_inte$Compounds)
length(a_1$Compounds)
length(a_3$Compounds)


# 纤维（f）
f1 <- Diff_1[Diff_1$panel%in%c("PvsW"),]%>%.[.$group=="Brownrice",]
f3 <- Diff_3[Diff_3$panel%in%c("PvsW"),]%>%.[.$group=="Brownrice",]
length(f1$Compounds)
length(f3$Compounds)
a
f_inte <- name[name$Compounds%in%intersect(f1$Compounds,f3$Compounds),]
f_all <- name[name$Compounds%in%c(f1$Compounds,f3$Compounds),]
f_1 <- name[name$Compounds%in%setdiff(f1$Compounds,f3$Compounds),]
f_3 <- name[name$Compounds%in%setdiff(f3$Compounds,f1$Compounds),]
write.table(f_all,"f_all.txt", sep = '\t', col.names = NA, quote = FALSE)
length(f_all$Compounds)
length(f_inte$Compounds)
length(f_1$Compounds)
length(f_3$Compounds)

# WT
wt1 <- Diff_1[Diff_1$panel%in%c("CvsWT"),]%>%.[.$group=="WildType",]
wt3 <- Diff_3[Diff_3$panel%in%c("CvsWT"),]%>%.[.$group=="WildType",]
length(wt1$Compounds)
length(wt3$Compounds)
wt_inte <- name[name$Compounds%in%intersect(wt1$Compounds,wt3$Compounds),]
wt_all <- name[name$Compounds%in%c(wt1$Compounds,wt3$Compounds),]
wt_1 <- name[name$Compounds%in%setdiff(wt1$Compounds,wt3$Compounds),]
wt_3 <- name[name$Compounds%in%setdiff(wt3$Compounds,wt1$Compounds),]
write.table(wt_all,"wt_all.txt", sep = '\t', col.names = NA, quote = FALSE)
length(wt_all$Compounds)
length(wt_inte$Compounds)
length(wt_1$Compounds)
length(wt_3$Compounds)
# metabolites origin
#metab_origin
Metabolites_Origin_br <- read.csv("tu/br_all/02_Origin_Analysis/Metabolites_Origin.csv")
Metabolites_Origin_br$Origin <- gsub("Drug related", "Others",Metabolites_Origin_br$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin_br$Origin <- factor(Metabolites_Origin_br$Origin,levels = c("Host","Microbiota","Co-Metabolism","Food related","Others"  ))
Metabolites_Origin_br$Class <- "B"
Metabolites_Origin_br$value <- 1
table(Metabolites_Origin_br$Origin)
Metabolites_Origin_wr <- read.csv("tu/wr_all/02_Origin_Analysis/Metabolites_Origin.csv")
Metabolites_Origin_wr$Origin <- gsub("Drug related", "Others",Metabolites_Origin_wr$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin_wr$Origin <- factor(Metabolites_Origin_wr$Origin,levels = c("Host","Microbiota","Co-Metabolism","Food related","Others"  ))
Metabolites_Origin_wr$Class <- "W"
Metabolites_Origin_wr$value <- 1
Metabolites_Origin_pr <- read.csv("tu/pr_all/02_Origin_Analysis/Metabolites_Origin.csv")
Metabolites_Origin_pr$Origin <- gsub("Drug related", "Others",Metabolites_Origin_pr$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin_pr$Origin <- factor(Metabolites_Origin_pr$Origin,levels = c("Host","Microbiota","Co-Metabolism","Food related","Others"  ))
Metabolites_Origin_pr$Class <- "P"
Metabolites_Origin_pr$value <- 1
Metabolites_Origin_s <- read.csv("tu/s_all/02_Origin_Analysis/Metabolites_Origin.csv")
Metabolites_Origin_s$Origin <- gsub("Drug related", "Others",Metabolites_Origin_s$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin_s$Origin <- factor(Metabolites_Origin_s$Origin,levels = c("Host","Microbiota","Co-Metabolism","Food related","Others"  ))
Metabolites_Origin_s$Class <- "S"
Metabolites_Origin_s$value <- 1
Metabolites_Origin_f <- read.csv("tu/f_all//02_Origin_Analysis/Metabolites_Origin.csv")
Metabolites_Origin_f$Origin <- gsub("Drug related", "Others",Metabolites_Origin_f$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin_f$Origin <- factor(Metabolites_Origin_f$Origin,levels = c("Host","Microbiota","Co-Metabolism","Food related","Others"  ))
Metabolites_Origin_f$Class <- "F"
Metabolites_Origin_f$value <- 1
Metabolites_Origin_a <- read.csv("tu/a_all/02_Origin_Analysis/Metabolites_Origin.csv")
Metabolites_Origin_a$Origin <- gsub("Drug related", "Others",Metabolites_Origin_a$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin_a$Origin <- factor(Metabolites_Origin_a$Origin,levels = c("Host","Microbiota","Co-Metabolism","Food related","Others"  ))
Metabolites_Origin_a$Class <- "A"
Metabolites_Origin_a$value <- 1

Metabolites_Origin_wt <- read.csv("tu/wt_all/02_Origin_Analysis/Metabolites_Origin.csv")
Metabolites_Origin_wt$Origin <- gsub("Drug related", "Others",Metabolites_Origin_wt$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin_wt$Origin <- factor(Metabolites_Origin_wt$Origin,levels = c("Host","Microbiota","Co-Metabolism","Food related","Others"  ))
Metabolites_Origin_wt$Class <- "WT"
Metabolites_Origin_wt$value <- 1
data <- rbind(Metabolites_Origin_br,Metabolites_Origin_wr)%>%
  rbind(.,Metabolites_Origin_pr)%>%
  rbind(.,Metabolites_Origin_s)%>%
  rbind(.,Metabolites_Origin_f)%>%
  rbind(.,Metabolites_Origin_a)%>%
  rbind(.,Metabolites_Origin_wt)
head(data)
data <- data%>%group_by(Class)%>%summarise(n=table(Origin))%>%as.data.frame()
data$Origin <- rep(c("Host", "Microbiota", "Co-Metabolism", "Food related", "Others"),7)
data$Origin <- factor(data$Origin,levels =c("Host", "Microbiota", "Co-Metabolism", "Food related", "Others") )
#View(data)
library(dplyr)
library(reshape2)
data$Class <- factor(data$Class,levels = c("S","B","W","P","F","A","WT"))

ggplot(data, aes(x=Class,y=n,fill=Origin)) +  
  geom_col(position = position_dodge(width = 0.8),width = 0.8) +  
  geom_text(aes(label=n),position = position_dodge(width = 0.8),vjust=-0.)+
  theme_classic()+scale_fill_manual(values = c("#A0B3A1","#4D7587","#85AEC9","#976276","#E1AB8D"))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+xlab("")+ylab("No. of Origin")

AOM_meta <- data[grep(c("Others"),data$Origin,invert = T),]%>%
  .[grep(c("Food related"),.$Origin,invert = T),]
write.table(AOM_meta,"AOM_meta.txt", sep = '\t', col.names = NA, quote = FALSE)



# bp
Co <- read.csv("tu/br_all/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
Micro <- read.csv("tu/br_all/03_Function_Analysis/MPEA_Results_Microbiota.csv")
H <- read.csv("tu/br_all/03_Function_Analysis/MPEA_Results_Host.csv")
Co$group <- "Co_Metabolism"
Micro$group <- "Micro_Metabolism"
H$group <- "Host_Metabolism"
bp <- rbind(Co,Micro)%>%rbind(.,H)
bp$logPvalue <- -log(bp$Pvalue)
bp <- bp[order(bp$Pvalue,decreasing = T),]
bp$Class <- "B"

# wp
Co <- read.csv("tu/wr_all/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
Micro <- read.csv("tu/wr_all/03_Function_Analysis/MPEA_Results_Microbiota.csv")
H <- read.csv("tu/wr_all/03_Function_Analysis/MPEA_Results_Host.csv")
Co$group <- "Co_Metabolism"
Micro$group <- "Micro_Metabolism"
H$group <- "Host_Metabolism"
wp <- rbind(Co,Micro)%>%rbind(.,H)
wp$logPvalue <- -log(wp$Pvalue)
wp <- wp[order(wp$Pvalue,decreasing = T),]
wp$Class <- "W"

# pp
Co <- read.csv("tu/pr_all/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
Micro <- read.csv("tu/pr_all/03_Function_Analysis/MPEA_Results_Microbiota.csv")
H <- read.csv("tu/pr_all/03_Function_Analysis/MPEA_Results_Host.csv")
Co$group <- "Co_Metabolism"
Micro$group <- "Micro_Metabolism"
H$group <- "Host_Metabolism"
pp <- rbind(Co,Micro)%>%rbind(.,H)
pp$logPvalue <- -log(pp$Pvalue)
pp <- pp[order(pp$Pvalue,decreasing = T),]
pp$Class <- "P"
#View(pp)
# sp
Co <- read.csv("tu/s_all/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
Micro <- read.csv("tu/s_all/03_Function_Analysis/MPEA_Results_Microbiota.csv")
H <- read.csv("tu/s_all/03_Function_Analysis/MPEA_Results_Host.csv")
Co$group <- "Co_Metabolism"
Micro$group <- "Micro_Metabolism"
H$group <- "Host_Metabolism"
sp <- rbind(Co,Micro)%>%rbind(.,H)
sp$logPvalue <- -log(sp$Pvalue)
sp <- sp[order(sp$Pvalue,decreasing = F),]
sp$Class <- "S"

# fp
Co <- read.csv("tu/f_all//03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
Micro <- read.csv("tu/f_all/03_Function_Analysis/MPEA_Results_Microbiota.csv")
H <- read.csv("tu/f_all/03_Function_Analysis/MPEA_Results_Host.csv")
Co$group <- "Co_Metabolism"
Micro$group <- "Micro_Metabolism"
H$group <- "Host_Metabolism"
fp <- rbind(Co,Micro)%>%rbind(.,H)
fp$logPvalue <- -log(fp$Pvalue)
fp <- fp[order(fp$Pvalue,decreasing = T),]
fp$Class <- "F"

# ap
Co <- read.csv("tu/a_all/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
Micro <- read.csv("tu/a_all/03_Function_Analysis/MPEA_Results_Microbiota.csv")
H <- read.csv("tu/a_all/03_Function_Analysis/MPEA_Results_Host.csv")
Co$group <- "Co_Metabolism"
Micro$group <- "Micro_Metabolism"
H$group <- "Host_Metabolism"
ap <- rbind(Co,Micro)%>%rbind(.,H)
ap$logPvalue <- -log(ap$Pvalue)
ap <- ap[order(ap$Pvalue,decreasing = T),]
ap$Class <- "A"

# wtp
Co <- read.csv("tu/wt_all/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
Micro <- read.csv("tu/wt_all/03_Function_Analysis/MPEA_Results_Microbiota.csv")
H <- read.csv("tu/wt_all/03_Function_Analysis/MPEA_Results_Host.csv")
Co$group <- "Co_Metabolism"
Micro$group <- "Micro_Metabolism"
H$group <- "Host_Metabolism"
wtp <- rbind(Co,Micro)%>%rbind(.,H)
wtp$logPvalue <- -log(wtp$Pvalue)
wtp <- wtp[order(wtp$Pvalue,decreasing = T),]
wtp$Class <- "WT"


data <- rbind(bp,wp)%>%rbind(.,pp)%>%rbind(.,sp)%>%rbind(.,fp)%>%rbind(.,ap)%>%rbind(.,wtp)
head(data)
data <- subset(data,logPvalue>3)
data$group <- factor(data$group,levels=c("Host_Metabolism","Micro_Metabolism","Co_Metabolism"))
data$Class <- factor(data$Class,levels = c("S","B","W","P","F","A","WT"))
#View(bp)
ggplot(data[data$Class=="WT",], aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#A0B3A1","#4D7587","#85AEC9"))+theme_classic()+
  ylab("")+xlab("log Pvalue")+scale_y_discrete(limits=unique(data[data$Class=="WT",]$Name))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 18,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))

ggplot(data[data$Class=="S",], aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#A0B3A1","#4D7587","#85AEC9"))+theme_classic()+
  ylab("")+xlab("log Pvalue")+scale_y_discrete(limits=rev(unique(data[data$Class=="S",]$Name)))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 18,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))

ggplot(data[data$Class=="B",], aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#A0B3A1","#4D7587","#85AEC9"))+theme_classic()+
  ylab("")+xlab("log Pvalue")+scale_y_discrete(limits=unique(data[data$Class=="B",]$Name))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 18,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))
ggplot(data[data$Class=="W",], aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#A0B3A1","#4D7587","#85AEC9"))+theme_classic()+
  ylab("")+xlab("log Pvalue")+scale_y_discrete(limits=unique(data[data$Class=="W",]$Name))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 18,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))

ggplot(data[data$Class=="P",], aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#A0B3A1","#4D7587","#85AEC9"))+theme_classic()+
  ylab("")+xlab("log Pvalue")+scale_y_discrete(limits=unique(data[data$Class=="P",]$Name))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 18,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))

ggplot(data[data$Class=="F",], aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#A0B3A1","#4D7587","#85AEC9"))+theme_classic()+
  ylab("")+xlab("log Pvalue")+scale_y_discrete(limits=unique(data[data$Class=="F",]$Name))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 18,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))

ggplot(data[data$Class=="A",], aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#A0B3A1","#4D7587","#85AEC9"))+theme_classic()+
  ylab("")+xlab("log Pvalue")+scale_y_discrete(limits=unique(data[data$Class=="A",]$Name))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 18,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))

colnames(logc)
logc$gene <- NULL
phData <- apply(logc[unique(c(Diff_1$Compounds,Diff_3$Compounds)),],1,function(x){tapply(x,sampleInfo$group,median)}) %>% t() %>% as.data.frame()
colnames(phData)
#phData <- phData[,c(1,6,2,7,3,8,4,9,5,10)]
n <- t(scale(t(phData)))
n[n>2]=2
n[n<-2]=-2
library(ComplexHeatmap)
library(circlize)

top_anno = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#83BFE3","#125AAF","#FFDD93","#F5C31E","#0ECBC5",
                                          "#0293A1","#E56A54","#F73E00","#D5EFA7","#56D75E")),
                       labels =c("WT_1","WT_3","C_1","C_3","B_1","B_3","W_1","W_3","P_1","P_3"),
                       labels_gp = gpar(col = "white", fontsize = 14)))
#left_anno= rowAnnotation(Panel = anno_block(gp = gpar(fill = c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F")),
#                                                 labels = levels(Diff$panel), 
#                                                 labels_gp = gpar(col = "white", fontsize = 12)))

select <- phData[rownames(phData)%in%c("LPC(22:5/0:0)","LPC(O-18:1/0:0)","LPC(22:6/0:0)",
                                       "Taurocholic acid","LPC(22:5/0:0)","Maritimetin",
                                       "Ergosine","Veratridine","Acetylvalerenolic acid","Methoprene acid",
                                       "Isosinomenine A","Epiafzelechin","Stylopine","Ajmalicine",
                                       "Deacetylvindoline","Mitragynine","Ursodeoxycholic Acid","L-Rhamnose",
                                       "Protocatechuic Aldehyde","4-Methylcatechol","butyrate","LPC(16:0/0:0)",
                                       "Phosphocholine","D-Turanose","D-(+)-sucrose","LPC(O-18:1/0:0)",
                                       "Deoxycholic acid","Glycochenodeoxycholic Acid","Dehydrolithocholic acid",
                                       "Proscillaridin A","Dehydrolithocholic acid","3,7-Di-O-methylquercetin",
                                       "4-Methylcatechol","D-Arabitol","Protocatechuic Aldehyde",
                                       "FFA(18:2)","FFA(20:0)","Indole-3-lactic acid"),]

select <- phData[rownames(phData)%in%c("LPC(16:0/0:0)",
                                       "Phosphocholine","D-Turanose","D-(+)-sucrose","LPC(O-18:1/0:0)",
                                       "Deoxycholic acid","Glycochenodeoxycholic Acid","Dehydrolithocholic acid",
                                       "Proscillaridin A","Dehydrolithocholic acid","3,7-Di-O-methylquercetin",
                                       "4-Methylcatechol","D-Arabitol","Protocatechuic Aldehyde",
                                       "FFA(18:2)","FFA(20:0)","Indole-3-lactic acid"),]
select <- phData[rownames(phData)%in%c("Indole-3-lactic acid"),]

dim(select)
gene_pos <- which(rownames(phData)  %in% rownames(select))
#gene_pos <- which(rownames(phData) %in% unique(rownames(select)))
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,
                                                 labels = rownames(select)))
col_fun = colorRamp2(c(-1, 0, 1), c("cornflowerblue","white","red"))
n=t(scale(t(phData)))
p1 <- Heatmap(n,col=col_fun,border = "black",
              show_row_names = FALSE,show_column_names = F,cluster_columns = FALSE,
              cluster_rows =T,
              column_split = factor(colnames(phData),levels=colnames(phData)),
              #              row_split = Diff$panel,
              top_annotation = top_anno,
              #              left_annotation = left_anno,
              right_annotation = row_anno,row_title = NULL,column_title = NULL,
              heatmap_legend_param = list(
                title = "scale",
                title_position = "leftcenter-rot"))
p1

head(Diff_1)
Diff_1$stage <- "1"
Diff_3$stage <- "3"
Diff <- rbind(Diff_1,Diff_3)
logc$gene <- rownames(logc)
data=melt(logc)
data$group <- paste0(strsplit2(data$variable,"_")[,1],"_",strsplit2(data$variable,"_")[,3])
data$group <- factor(data$group,levels =c("WT_1","WT_3","C_1","C_3","B_1","B_3","W_1","W_3","P_1","P_3") )
Data_summary <- summarySE(data, measurevar="value", groupvars=c("group","variable"))
View(ggdata)
ggdata <- data[data$gene=="FFA(18:2)",]
ggplot(ggdata,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#83BFE3","#125AAF","#FFDD93","#F5C31E","#0ECBC5",
                             "#0293A1","#E56A54","#F73E00","#D5EFA7","#56D75E"))+
  geom_signif(comparisons = list(c("WT_1", "C_1"),c("WT_3", "C_3"),
                                 c("C_1","B_1"),c("C_3","B_3"),
                                 c("C_1","W_1"),c("C_3","W_3"),
                                 c("C_1","P_1"),c("C_3","P_3")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative concentration of L-Tryptophan")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

ggdata <- data[data$gene=="AA",]
ggplot(ggdata,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#83BFE3","#125AAF","#FFDD93","#F5C31E","#0ECBC5",
                             "#0293A1","#E56A54","#F73E00","#D5EFA7","#56D75E"))+
  geom_signif(comparisons = list(c("WT_1", "C_1"),c("WT_3", "C_3"),
                                 c("C_1","B_1"),c("C_3","B_3"),
                                 c("C_1","W_1"),c("C_3","W_3"),
                                 c("C_1","P_1"),c("C_3","P_3")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative concentration of L-Histidine")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

## 维B  Folic-acid
ggdata <- data[data$gene=="L-kynurenine",]
ggplot(ggdata,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#83BFE3","#125AAF","#FFDD93","#F5C31E","#0ECBC5",
                             "#0293A1","#E56A54","#F73E00","#D5EFA7","#56D75E"))+
  geom_signif(comparisons = list(c("WT_1", "C_1"),c("WT_3", "C_3"),
                                 c("C_1","B_1"),c("C_3","B_3"),
                                 c("C_1","W_1"),c("C_3","W_3"),
                                 c("C_1","P_1"),c("C_3","P_3")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative concentration of 
Tryptamine")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

################################# Correlation #####################################
AOM_meta <- read.delim("AOM_meta.txt")
dim(AOM_meta)
B <- AOM_meta[AOM_meta$Class=="B",]
head(B)
intersect(B$Name,wt_all$Compounds)
health <- intersect(AOM_meta$Name,wt_all$Compounds)

a <- setdiff(AOM_meta$Name,name$Compounds)
AOM_meta[AOM_meta$HMDBID%in%a,]
head(phData)
mean_metab <- phData
write.table(mean_metab,"mean_metab.txt", sep = '\t', col.names = NA, quote = FALSE)

diff_metab <- mean_metab[rownames(mean_metab)%in%Apc_meta$Name,]
mean_micro <- apply(otu_order,1,function(x){tapply(x,Black_riceInfo$group,mean)}) %>% t() %>% as.data.frame()
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
meta <- data.frame(Metabolite_sampleID=a,Microbiome_sampleID=a,
                   Grouping=paste0(strsplit2(a,"_")[,1],"_",strsplit2(a,"_")[,3]))
head(meta)
data1 <- diff_metab[,a]
head(name)
data1$Compounds <- rownames(data1)
data1 <- merge(data1,name,all.x=T,by="Compounds")
head(data1)
data2 <- diff_micro[,a]
write.table(meta,"meta.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(data1,"data_metab.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(data2,"data_micro.txt", sep = '\t', col.names = NA, quote = FALSE)

library(psych)
cor <-corr.test(t(diff_metab[,a]),t(diff_micro[,a]), method = "spearman",adjust= "none")
cor <-corr.test(t(diff_metab),t(diff_micro), method = "spearman",adjust= "none")

#cor <-corr.test(diff_microb, diff_metab, method = "pearson",adjust= "none")
## 提取相关性、p值
cor
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
a=unique(Apc_meta[Apc_meta$Class%in%c("B","W","P","A","F"),]$Name)
length(a)
dim(cmt)
a=intersect(a,rownames(cmt))
data <- cmt[a,]
library(ComplexHeatmap)
library(circlize)

left_anno= rowAnnotation(Panel = anno_block(gp = gpar(fill = c("#C2D6E4","#FBE195","#CFD9D0")),
                                            labels = c("B","W","P"), 
                                            labels_gp = gpar(col = "white", fontsize = 12)))
select <- data[rownames(data)%in%c(health,"L-Histidine","L-Tryptophan","Indole-3-lactic acid",
                                   "Estriol","Biopterin","D-Tagatose","4-Hydroxybenzoic Acid",
                                   "14,15-DHET","Riboflavin","Phosphocholine","Malonicacid"),]
select
#"L-Alanine","L-Glutamic Acid","D-Galactose","Myoinositol",
#"ADP-ribose","Adenosine","Hypoxanthine","Xanthine"
select <- select[!duplicated(rownames(select)),]
dim(select)
gene_pos <- which(rownames(data)  %in% rownames(select))
gene_pos <- which(unique(rownames(data))%in% unique(rownames(select)))
col_fun = colorRamp2(c(-1,0,1), c("cornflowerblue","white","red"))
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,
                                                 labels = unique(rownames(select))))
p1 <- Heatmap(data[,c(1,2,7,9,10,3,4,5,6,8)],col=col_fun,border = "black",
              show_row_names = FALSE,show_column_names =T,
              cluster_columns = F,cluster_rows = T,
              #left_annotation = left_anno,
              right_annotation = row_anno,
              row_title = NULL,column_title = NULL,
              heatmap_legend_param = list(title = "scale",
                                          title_position = "leftcenter-rot"))
p1


pheatmap::pheatmap(cmt,fontsize = 20,show_colnames = T,legend_breaks = c(-1,-0.5,0,0.5,1),display_numbers = t(pmt),
                   color=mycol)





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

