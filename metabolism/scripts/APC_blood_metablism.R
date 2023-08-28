setwd("F:/A_MetaboAnalysis/Black_rice/Black_rice_metab/blood_metab/")
#################################### Apc BR #########################################
library(readxl)
library(tidyverse)
library(limma)
ALL_sample_data <- read_delim("ALL_data.txt")
View(ALL_sample_data)
colnames(ALL_sample_data)
QC <- ALL_sample_data[,26:28]
head(QC)
ALL_meta <- data.frame(Compounds=ALL_sample_data$Compounds,ALL_sample_data[,12:25])
colnames(ALL_meta)
ALL_meta[1:5,1:5]
rownames(ALL_meta)<- ALL_meta$Compounds
ALL_meta$Compounds <- NULL
head(ALL_meta)
sampleInfo <- data.frame(substr(colnames(ALL_meta),1,1))
colnames(sampleInfo) <- c("group")
rownames(sampleInfo) <- colnames(ALL_meta)
sampleInfo$group <- factor(sampleInfo$group,levels = c("C","B"))
dim(sampleInfo)
colnames(ALL_meta)
logc <- log(ALL_meta+1)
head(logc)
rownames(ALL_meta)[1:5]
dim(logc)
scale_logc <- t(scale(t(logc)))
library(PCAtools)
pca <- pca(logc, metadata =sampleInfo)
biplot(pca, x = 'PC1', y = 'PC2') #可以看到这两个成分对样品的解释度
#screen(p) #可以看到所有成分对样品的解释
##将pca与关联样本信息表 
pca$variance
pca_rotated_plus <- cbind(pca$rotated,sampleInfo)
library(ggsci)
pca_rotated_plus$group <- factor(pca_rotated_plus$group,levels=c("C","B"))

ggplot(pca_rotated_plus,aes(x = PC1 , y = PC2))+
  geom_point(size = 8,aes(#shape = group, 
    fill = group,color = group)) +
  stat_ellipse(aes(color = group,fill=group),linetype = 'dashed',size = 1, show.legend = FALSE) + #添加分组椭圆
  labs(x = 'PC1 (23.72%)',y = 'PC2 (14.89%)') + 
  scale_fill_manual(values = c("#D14424","#0094A5"))+
  scale_color_manual(values = c("#D14424","#0094A5"))+
  theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.position = "top",legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))

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
View(logc)

Pvalue<-c(rep(0,nrow(logc)))
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,1:7])==0&&sd(logc[i,8:14])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,1:7]),as.numeric(logc[i,8:14]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(ALL_meta[i,1:7]))/mean(as.numeric(ALL_meta[i,8:14])) 
  }}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)

# VIP value 
library(ropls) 
colnames(scale_logc)
oplsda <-  opls(as.data.frame(t(scale_logc)), 
                sampleInfo$group, orthoI = 1) 
vip<- getVipVn(oplsda)%>%as.data.frame(.)
head(vip)

#提取样本在 OPLS-DA 轴上的位置
sample.score = oplsda@scoreMN %>%  #得分矩阵
  as.data.frame() %>%
  mutate(p1=p1,o1 = oplsda@orthoScoreMN[,1]) #正交矩阵

head(sample.score)#查看
sample.score$group <- substr(rownames(sample.score),1,1)

ggplot(sample.score, aes(p1, o1, color = group)) +
  geom_point(size = 8,aes(#shape = group, 
    fill = group,color = group)) +
  stat_ellipse(aes(color = group,fill=group),linetype = 'dashed',size = 1, show.legend = FALSE) + #添加分组椭圆
  labs(x = 'P1(25.5%)',y = 'to1') + 
  scale_fill_manual(values = c("#D14424","#0094A5"))+
  scale_color_manual(values = c("#D14424","#0094A5"))+
  theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.position = "top",legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))


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
write.table(res_MF,"res_MF.txt", sep = '\t', col.names = NA, quote = FALSE)

meta <- meta[meta$Metabolite%in%c("Sinapic acid","Levan","Indole","ChlorogenicAcid",
                                  "Biliverdin","DTMP","Folinic acid"),]
head(res_MF)
dim(meta)
table(res_MF$type)
p=ggplot(res_MF,aes(x = log2FC,y = -log10(FDR))) +
  geom_point(aes(color = type,size=VIP),alpha = 0.5,size = 3) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1,
        # 标题居中
        plot.title = element_text(),
        axis.text = element_text(color = 'black',size=24),axis.title.x  = element_text(color = 'black',size=24),
        axis.title = element_text(color = 'black'),axis.title.y   = element_text(color = 'black',size=24),
        legend.text = element_text(size=14),legend.position = "right") +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#DA1212','noSig'='grey','down'='#3E7C17'),
                     # legend labels
                     label = c('up'='up (C, num=113)','noSig'='noSig (num=764)','down'='down (B, num=110)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('Diet CvsB')
p
p +	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))





mix_select <- mix_select[order(mix_select$log2FC),]
mix_select$group <- 0
mix_select$group[which(mix_select$log2FC>0)] <- "Control"
mix_select$group[which(mix_select$log2FC<0)] <- "Blackrice"
Diff_C_B <- mix_select
Diff_C_B$Compounds <- rownames(Diff_C_B)
dim(Diff_C_B)
write.table(Diff_C_B,"Diff_C_B.txt", sep = '\t', col.names = NA, quote = FALSE)
View(ALL_sample_data[ALL_sample_data$Compounds%in%Diff_C_B[Diff_C_B$group=="Blackrice",]$Compounds,])
write.table(ALL_sample_data[ALL_sample_data$Compounds%in%Diff_C_B[Diff_C_B$group=="Blackrice",]$Compounds,][,c(2,30)],"metab_name.txt", sep = '\t', col.names = NA, quote = FALSE)
colnames(Diff_C_B)
Mean <- apply(Diff_C_B[,1:14],1,function(x){tapply(x,substr(colnames(Diff_C_B[,1:14]),1,1),mean)}) %>% t() %>% as.data.frame()
head(Mean)




Metabolites_Origin<- read.csv("MetOrigin_20221110224441/02_Origin_Analysis/Metabolites_Origin.csv")
Metabolites_Origin$Origin <- gsub("Drug related", "Others",Metabolites_Origin$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin$Origin <- factor(Metabolites_Origin$Origin,levels = c("Host","Microbiota","Co-Metabolism","Food related","Others"  ))
Metabolites_Origin$Class <- "CvsB"
Metabolites_Origin$value <- 1
data <- as.data.frame(table(Metabolites_Origin$Origin))
colnames(data) <- c("Origin","n")
ggplot(data, aes(x=Origin,y=n,fill=Origin)) +  
  geom_col(position = position_dodge(width = 0.8),width = 0.5) +  
  geom_text(aes(label=n),position = position_dodge(width = 0.5),vjust=-0.)+
  theme_classic()+scale_fill_manual(values = c("#A0B3A1","#4D7587","#85AEC9","#976276","#E1AB8D"))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold",angle = 45,hjust = 1),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+xlab("")+ylab("N of Origin")
library(RColorBrewer)
library(dplyr)
library(graphics)
library(ggplot2)
df <-arrange(data,n)
labs <- paste0(df$Origin," \n(", round(df$n/sum(df$n)*100,2), "%)") #标签
lab <- paste0(round(df$n/sum(df$value)*100,2), "%") #标签
pie(df$n,labels=labs, init.angle=90,col =  c("#A0B3A1","#4D7587","#85AEC9","#976276","#E1AB8D"),
    border="black")
Co <- read.csv("MetOrigin_20221110224441/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
Micro <- read.csv("MetOrigin_20221110224441/03_Function_Analysis/MPEA_Results_Microbiota.csv")
H <- read.csv("MetOrigin_20221110224441/03_Function_Analysis/MPEA_Results_Host.csv")
Co$group <- "Co_Metabolism"
Micro$group <- "Micro_Metabolism"
H$group <- "Host_Metabolism"
data <- rbind(Co,Micro)%>%rbind(.,H)
data$logPvalue <- -log(data$Pvalue)
data <- data[order(data$Pvalue,decreasing = T),]

data <- subset(data,logPvalue>2.6)
data$group <- factor(data$group,levels=c("Host_Metabolism","Micro_Metabolism","Co_Metabolism"))

#View(bp)
ggplot(data, aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#A0B3A1","#4D7587","#85AEC9"))+theme_classic()+
  ylab("")+xlab("-log Pvalue")+scale_y_discrete(limits=unique(data$Name))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 18,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))
