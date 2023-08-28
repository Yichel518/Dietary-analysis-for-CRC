setwd("F:/A_MetaboAnalysis/Black_rice/")
#devtools::install_local("F:/A_aging_lnc/LDM-master/LDM-master/LDM_5.0.tar.gz")
library("LDM")
#remotes::install_github("https://github.com/xielab2017/EasyMicroPlot",subdir='Version_0.5')
library(EasyMicroPlot)
#install.packages("ggiraph")
#install.packages("rlang")
#install.packages("vctrs")
library(LDM)
library(reshape2)
library(dplyr)
# 种 species
Black_rice_species <- read.delim("F:/A_MetaboAnalysis/Black_rice/Apc_all_merge_abundance.txt",row.names = 1 ,header=T)
Black_rice_species <- Black_rice_species[,-c(53:112)]
rownames(Black_rice_species)
library(reshape2)
library(dplyr)
library(tidyverse)
library(limma)
colnames(Black_rice_species)
dim(Black_rice_species)

Black_riceInfo <- strsplit2(colnames(Black_rice_species),"_")%>%as.data.frame()
Black_riceInfo$V3 <- substr(Black_riceInfo$V3,1,3)
colnames(Black_riceInfo) <- c("condition","sex","stage")
rownames(Black_riceInfo) <- colnames(Black_rice_species)
Black_riceInfo$group <- paste0(Black_riceInfo$condition,"_",Black_riceInfo$stage)
dim(Black_rice_species)
Black_riceInfo$condition <- factor(Black_riceInfo$condition,levels = c("WT","C","B"))
Black_riceInfo$group <- factor(Black_riceInfo$group,
                               levels = c("WT_14w","WT_22w","C_14w","C_22w",
                                          "B_14w","B_22w"))
Black_riceInfo <- Black_riceInfo[order(Black_riceInfo$group),]
Black_rice_species <- Black_rice_species[,rownames(Black_riceInfo)]
colnames(Black_rice_species)
# 通常我们只关注高丰度且显著差异的，按每个OTU的总丰度排序
Black_rice_species$sum <- rowSums(Black_rice_species)
# 按每行的和降序排列
Black_rice_species_order <- Black_rice_species[order(Black_rice_species$sum, decreasing=TRUE), ]
otu <- Black_rice_species_order
otu$sum <- rowSums(otu)
otu_order <- otu[order(otu$sum, decreasing=TRUE), ]
colnames(otu)
otu_order$sum <- NULL
# Correlation spearman
pheatmap::pheatmap(cor(otu_order,method = "spearman"),cluster_rows = F,cluster_cols = F)
dim(otu_order[1:20,])
mat1 <- otu_order[1:20,]
head(mat1)
mat2 <- otu_order[1:30,]
dim(Black_riceInfo)
tmp1 = apply(mat1,1,function(x){tapply(x,Black_riceInfo$group,mean)}) %>% t() %>% as.data.frame()
tmp2 = apply(mat2,1,function(x){tapply(x,Black_riceInfo$group,median)}) %>% t() %>% as.data.frame()
tmp = apply(otu_order,1,function(x){tapply(x,Black_riceInfo$group,median)}) %>% t() %>% as.data.frame()
dim(tmp)
n=t(scale(t(tmp)))
n[n>2]=2
n[n<-2]=-2
set.seed(1)
others <-apply(tmp[!rownames(tmp)%in%rownames(tmp1),],2,sum)%>%as.data.frame(.)
colnames(others) <- "Others"
others <- as.data.frame(t(others))
mat <- rbind(tmp1,others)
mat
ggdata = mat %>% rownames_to_column(var = "species") %>% gather(key = "group",value = "value", -"species")
head(ggdata)
unique(ggdata$species)
ggdata$species <- factor(ggdata$species,levels = unique(ggdata$species))
head(ggdata)
ggdata$group<- factor(ggdata$group,levels=c("WT_14w","C_14w","B_14w","WT_22w","C_22w","B_22w"))
ggdata$condition<- factor(ggdata$group,levels=c("WT","C","B"))
unique(ggdata$species)
ggplot(ggdata) +				
  geom_bar(aes(x = group,y=value,fill = species),stat = "identity",position = "fill",colour= 'black') +						
  scale_fill_manual(values  = c("#6686c6","#a0c1db","gray","#E65F92","#c77364","#ce8f5c",						
                                "#7bac80","#75a5c0","#b5181a","#b72d9e",						
                                "#e4cccb","#f6f0eb","#e8c755","#d5d456","#cfe74f","#39cf40",						
                                "#3e6926","#0a0883","#49ddd0","#e0f8f7",						
                                "#651d91","#9d3eaf","#b9a3d0","#5b5a5b","#8f8f8f"))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+	
  xlab("")+ylab("% of species")+guides(fill = guide_legend( ncol = 1, byrow = TRUE))						


library(vegan)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(tidyr)
library(ggpubr)

#######################获得多种统计alpha β多样性统计量###########################
library(vegan)
library(picante)      

alpha_diversity <- function(x, tree = NULL) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  #保留四位小数
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_Coverage <- sprintf("%0.4f", goods_Coverage)
  result <- data.frame(observed_species, ACE,Chao1, 
                       Shannon, Simpson, goods_Coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    
    result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson,
                         PD_whole_tree ,goods_Coverage)
  }
  
  
  result
}

count<- read.delim("F:/A_MetaboAnalysis/Black_rice/merged_counts_table_species.txt", row.names=1)
colnames(count)
dim(count)
count <- count[,-c(30:89)]
meta <- strsplit2(colnames(count),"[.]")[,1]%>%strsplit2(.,"_")%>%as.data.frame()
colnames(meta) <- c("condition","sex","stage")
rownames(meta) <- colnames(count)
meta$group <- paste0(meta$condition,"_",meta$stage)
dim(Black_rice_species)
meta$condition <- factor(meta$condition,levels = c("WT","C","B"))
meta$group <- factor(meta$group,levels = c("WT_14w","C_14w","B_14w","WT_22w","C_22w","B_22w"))
meta <- meta[order(meta$group),]
count <- count[,rownames(meta)]

otu_t <- t(count)
rownames(otu_t)

alpha <- alpha_diversity(otu_t)
alpha$ACE<-NULL
alpha$goods_Coverage <- NULL

alpha$group <- meta$group
alpha$condition <- meta$condition
alpha$sex <- meta$sex
head(alpha)
otu_diversity_long <- melt(alpha,id.vars = c("group","condition","sex")) 
otu_diversity_long$group <- factor(otu_diversity_long$group,levels=c("WT_14w","C_14w","B_14w",
                                                                     "WT_22w","C_22w","B_22w"))
head(otu_diversity_long)
otu_diversity_long$value <- as.numeric(otu_diversity_long$value)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
head(otu_diversity_long)
Data_summary <- summarySE(otu_diversity_long, measurevar="value", groupvars=c("group","condition","variable"))

library(rstatix)
otu_diversity_long$group <- factor(otu_diversity_long$group,levels =c("WT_14w","C_14w","B_14w",
                                                                      "WT_22w","C_22w","B_22w"))

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}
#去除outlier
#library(dplyr)
data <- otu_diversity_long[otu_diversity_long$variable=="Simpson",]
data1 <- data %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()   

data2 <- data %>% 
  group_by(condition) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
otu_diversity_long$group
ggplot(otu_diversity_long[otu_diversity_long$variable=="Simpson",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(notch = F,lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Simpson diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

ggplot(otu_diversity_long[otu_diversity_long$variable=="Simpson",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(notch = F,lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+ylab("Simpson diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')



ggplot(data2,aes(x = condition,y = value,fill=condition))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(notch = F,lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Simpson diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')






data <- otu_diversity_long[otu_diversity_long$variable=="Shannon",]
data1 <- data %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup() 

data2 <- data %>% 
  group_by(condition) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup() 


ggplot(data1,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Shannon diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')


ggplot(data2,aes(x = condition,y = value,fill=condition))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(notch = F,lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#E29997","#C3B3D5","#EBBE70","#C1DF8C","#B4CEE3"))+
  geom_signif(comparisons = list(c("WT", "C"),c("C","B"),c("WT","B")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Shannon diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

data <- otu_diversity_long[otu_diversity_long$variable=="observed_species",]
data1 <- data %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup() 

data2 <- data %>% 
  group_by(condition) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup() 
ggplot(data1,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("No. observed OTUs")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

ggplot(data2,aes(x = condition,y = value,fill=condition))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(notch = F,lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#E29997","#C3B3D5","#EBBE70","#C1DF8C","#B4CEE3"))+
  geom_signif(comparisons = list(c("WT", "C"),c("C","B"),c("WT", "B")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("No. observed OTUs")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')



data <- otu_diversity_long[otu_diversity_long$variable=="Chao1",]
data1 <- data %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup() 

data2 <- data %>% 
  group_by(condition) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup() 

ggplot(otu_diversity_long[otu_diversity_long$variable=="Chao1",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Chao1 diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

ggplot(data2,aes(x = condition,y = value,fill=condition))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(notch = F,lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#E29997","#C3B3D5","#EBBE70","#C1DF8C","#B4CEE3"))+
  geom_signif(comparisons = list(c("WT", "C"),c("C","B"),c("C", "W"),c("C","P"),
                                 c("B","W"),c("W","P"),c("B", "P")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Chao1 diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')



data <- otu_order[rownames(otu_order)%in%c("Bacteroides_uniformis",
                                           "Akkermansia_muciniphila",
                                           "Lactobacillus_johnsonii",
                                           "Escherichia_coli",
                                           "Lactobacillus_lactis",
                                           "Lactobacillus_murinus","Lactobacillus_reuteri"),]
data$species <- rownames(data)
ggdata <- melt(data)
ggdata
meta <- strsplit2(ggdata$variable,"[.]")[,1]%>%strsplit2(.,"_")%>%as.data.frame()
ggdata$group <- paste0(meta[,1],"_",meta[,3])%>%
  factor(.,levels = c("WT_14w","C_14w","B_14w",
                      "WT_22w","C_22w","B_22w"))
data <- ggdata[ggdata$species=="Bacteroides_uniformis",]
data1 <- data %>% 
  group_by(group) %>% 
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
  ylab("Relative abundance")+ggtitle("B.uniformis")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))

library(dplyr)

data <- ggdata[ggdata$species=="Lactobacillus_reuteri",]
data1 <- data %>% 
  group_by(group) %>% 
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
  ylab("Relative abundance")+ggtitle("L.reuteri")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))


data <- ggdata[ggdata$species=="Akkermansia_muciniphila",]
data1 <- data %>% 
  group_by(group) %>% 
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
  ylab("Relative abundance")+ggtitle("A.muciniphila")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))



data <- ggdata[ggdata$species=="Lactobacillus_johnsonii",]
data1 <- data %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup() 
ggplot(na.omit(data),aes(x = group,y = value,fill=group))+
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
  ylab("Relative abundance")+ggtitle("L.johnsonii")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))



data <- ggdata[ggdata$species=="Escherichia_coli",]
data1 <- data %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup() 
ggplot(na.omit(data),aes(x = group,y = value,fill=group))+
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
  ylab("Relative abundance")+ggtitle("E.coli")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))


####################################beita#########################################
?pcoa
library(ape)
metaphlan_anno = read.table("F:/A_MetaboAnalysis/metaphlan_anno.txt",header = T) 
metaphlan_tree = ape::read.tree("F:/A_MetaboAnalysis/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk")
metaphlan_tree$tip.label <- metaphlan_tree$tip.label %>% gsub(".+\\|s__","",.)
my_tree = ape::keep.tip(metaphlan_tree,intersect(rownames(count),metaphlan_tree$tip.label))
#BiocManager::install("ggdendro")
library(ggplot2)
library(ggdendro)
library(phyloseq)
library(tidyverse)
library(ggfun)
#BiocManager::install("amplicon")
library(ggtree)
library(ggstance)
library(phyloseq)
library(vegan)
dim(sampleInfo)
#meta <- strsplit2(colnames(count),"[.]")[,1]%>%strsplit2(.,"_")%>%as.data.frame()
#colnames(meta) <- c("condition","sex","stage")
#rownames(meta) <- colnames(count)
#meta$group <- paste0(meta$condition,"_",meta$stage)

#meta$condition <- factor(meta$condition,levels = c("WT","C","B"))
#meta$group <- factor(meta$group,
#                     levels = c("WT_14w","C_14w","B_14w",
#                                "WT_22w","C_22w","B_22w"))
#
#meta <- meta[order(meta$group),]
#count <- count[,rownames(meta)]
my_phy = phyloseq(otu_table(as.matrix(count),taxa_are_rows = TRUE),
                  sample_data(meta),phy_tree(my_tree),tax_table(metaphlan_anno %>% as.matrix()))  #这里是否要as.matrix很严格

weight_unifrac= phyloseq::distance(my_phy,method = "wunifrac") %>% as.matrix()
unweight_unifrac = phyloseq::distance(my_phy,method = "unifrac") %>% as.matrix()
bray = phyloseq::distance(my_phy,method = "bray") %>% as.matrix() %>% as.data.frame()

#合并展示
res = pcoa(bray,correction = "cailliez",rn = NULL)
biplot.pcoa(res)
res$vectors
pcoachoose = c("Axis.1","Axis.2")
ggdata = res$vectors[,pcoachoose] %>%as.data.frame()
dim(ggdata)
ggdata$group <- factor(meta$group,levels = c("WT_14w","C_14w","B_14w",
                                             "WT_22w","C_22w","B_22w"))

ggplot(data = ggdata,aes(x=get(pcoachoose[1]),y=get(pcoachoose[2]),fill=group,color = group))+						
  geom_point(size=8) +
  stat_ellipse(type = "t", linetype = 1, show.legend = FALSE,aes(color = group))+						
  ylab("Bray-Curtis diversity") +theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+						
  scale_color_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+						
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+						
  xlab("PCoA1 (27.8%)")+ylab("PCoA2 (13.7%)")


dune.dist <- vegdist(t(count), method="bray", binary=F)
dune.div <- adonis2(dune.dist ~ group, data = meta, permutations = 999, method="bray")
dune.div
dispersion <- betadisper(dune.dist, group=meta$group)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse


library(pairwiseAdonis)

dune.pairwise.adonis <- pairwise.adonis(x=bray, 
                                        factors=meta$group, 
                                        sim.function = "vegdist",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)
dune.pairwise.adonis

dune_pcoa <- cmdscale(bray, k=3, eig=T)
dune_pcoa_points <- as.data.frame(dune_pcoa$points)
sum_eig <- sum(dune_pcoa$eig)
eig_percent <- round(dune_pcoa$eig/sum_eig*100,1)
eig_percent
set.seed(1)

sub_dist <- list()
for (group in unique(meta$condition)) { 
  row_group <- rownames(meta[which(meta$condition==group),])
  sample_group <- t(count)[row_group,]
  sub_dist[[group]] <- bray[ row_group, row_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups<- melt(sub_dist)
df.bray <- braygroups[complete.cases(braygroups), ]
df.bray$L1 <- factor(df.bray$L1, levels=c("WT","C","B","W","P"))

df.bray$L1 <- factor(df.bray$L1, levels=c("WT_14w","C_14w","B_14w",
                                          "WT_22w","C_22w","B_22w"))

head(df.bray)

ggplot(df.bray, aes(x=L1, y=value, fill=L1)) +
  geom_jitter() +
  scale_colour_ordinal()+
  # geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=0.6,position=position_dodge(0.9))+
  theme(legend.position="none") +
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+						
  ylab("Bray-Curtis diversity") +theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+xlab("")+
  geom_signif(comparisons = list(c("WT","C"),c("C","B"),c("WT","B")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)

sub_dist <- list()
for (group in unique(meta$group)) { 
  row_group <- rownames(meta[which(meta$group==group),])
  sample_group <- t(count)[row_group,]
  sub_dist[[group]] <- bray[ row_group, row_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups<- melt(sub_dist)
df.bray <- braygroups[complete.cases(braygroups), ]
df.bray$L1 <- factor(df.bray$L1, levels=c("WT_14w","C_14w","B_14w",
                                          "WT_22w","C_22w","B_22w"))
data <- df.bray %>% 
  group_by(L1) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()  

ggplot(df.bray, aes(x=L1, y=value, fill=L1)) +
  geom_jitter() +
  scale_colour_ordinal()+
  # geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=0.6,position=position_dodge(0.9))+
  theme(legend.position="none")+					
  ylab("Bray-Curtis diversity") +theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45,hjust=1,face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+xlab("")+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)



######## 可视化 #############
library(ape)          #做PCoA
library(RColorBrewer) #配色
#BiocManager::install("usedist")
library(usedist)      #PERMANOVA test时处理距离矩阵

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
# 读取表达矩阵


ggdata$Axis.1
dim(ggdata)
data1 <- ggdata %>% 
  group_by(group) %>% 
  mutate(Axis.1 = remove_outliers(Axis.1),Axis.2 = remove_outliers(Axis.2)) %>% 
  ungroup()
dim(data1)
plotinfo <- na.omit(data1)
plotinfo$stage <- factor(strsplit2(plotinfo$group,"_")[,2],levels = c("14w","22w"))
plotinfo$group <- factor(plotinfo$group,levels = c("WT_14w","C_14w","B_14w",
                                                   "WT_22w","C_22w","B_22w"))
table(plotinfo$group)
plotinfo2 <- NULL
tmp <- plotinfo[plotinfo$group == "WT_14w",]
tmp
for (i in unique(plotinfo$group)) {
  tmp <- plotinfo[plotinfo$group == i,] # 取出当前亚型当前来源下的数据
  avgx <- mean(tmp$Axis.1) # 计算横坐标均值
  avgy <- mean(tmp$Axis.2) # 计算纵坐标均值
  sdx <- sd(tmp$Axis.1) # 计算横坐标标准差
  sdy <- sd(tmp$Axis.2) # 计算纵坐标标准差
  
  plotinfo2 <- rbind.data.frame(plotinfo2,
                                data.frame(group = i,
                                           shape = unique(ifelse(unique(tmp$stage) == "14w", # 如果是TCGA就实心圆
                                                                 "opened","opened")), # 如果是Yau就为空心圆
                                           label = unique(tmp$group),
                                           avgx = avgx, # 添加圆的x位置
                                           avgy = avgy, # 添加圆的y位置
                                           sdx = sdx, # 添加圆的水平标准差
                                           sdy = sdy, # 添加圆的垂直标准差
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}




plotinfo2$color <- c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00")

# 绘图
?par
par(las = 0)
par(bty="l", mgp = c(2.5,.33,0), mar=c(5,5,5,5)+.1, las=1, tcl=-.25)
# 根据上面的基础图像，调整横轴纵轴的宽度，生成一张空白图像
?plot
plot(NULL,NULL,
     xlab = "PcoA:27.8%",
     ylab = "PcoA:13.7%",
     xlim = c(-0.35,0.3),
     ylim = c(-0.25,0.4),cex.axis=1.5,cex.lab=2)
?plot
# 先产生误差线，这样可以被后面的圆挡住
for (i in 1:nrow(plotinfo2)) {
  tmp <- plotinfo2[i,]
  
  # 产生横向误差线
  lines(x = c(tmp$avgx - tmp$sdx,
              tmp$avgx + tmp$sdx),
        y = c(tmp$avgy, tmp$avgy),
        lty = ifelse(tmp$shape == "closed",1,2), # 如果是closed就是实线，否则为虚线
        col = tmp$color,
        lwd = 2)
  
  # 产生横向误差线的封口线
  lines(x = c(tmp$avgx - tmp$sdx,
              tmp$avgx - tmp$sdx),
        y = c(tmp$avgy - 0.01, # 宽度根据情况调整
              tmp$avgy + 0.01),
        col = tmp$color,
        lwd = 2)
  lines(x = c(tmp$avgx + tmp$sdx,
              tmp$avgx + tmp$sdx),
        y = c(tmp$avgy - 0.01,
              tmp$avgy + 0.01),
        col = tmp$color,
        lwd = 2)
  
  # 产生纵向误差线
  lines(x = c(tmp$avgx, tmp$avgx),
        y = c(tmp$avgy - tmp$sdy,
              tmp$avgy + tmp$sdy),
        lty = ifelse(tmp$shape == "closed",1,2),
        col = tmp$color,
        lwd = 2)
  
  # 产生纵向误差线的封口线
  lines(x = c(tmp$avgx -0.01,
              tmp$avgx + 0.01),
        y = c(tmp$avgy + tmp$sdy,
              tmp$avgy + tmp$sdy),
        col = tmp$color,
        lwd = 2)
  lines(x = c(tmp$avgx - 0.01,
              tmp$avgx + 0.01),
        y = c(tmp$avgy - tmp$sdy,
              tmp$avgy - tmp$sdy),
        col = tmp$color,
        lwd = 2)
}

# 后添加圆，以挡住误差线
points(plotinfo2$avgx,
       plotinfo2$avgy,
       pch = ifelse(plotinfo2$shape == "closed",19,21), # 如果为closed就是实心圆，否则为空心圆
       bg = ifelse(plotinfo2$shape == "closed",plotinfo2$color,"white"), # 填充背景色，如果为closed就是为该颜色，否则为白色
       col = plotinfo2$color, # 填充边框颜色
       lwd = 2, # 边框粗细
       cex = 8) # 圆的大小

# 添加文本
text(plotinfo2$avgx,
     plotinfo2$avgy,
     plotinfo2$label,
     col = ifelse(plotinfo2$shape == "closed","white",plotinfo2$color), # 如果是实心圆文字为白色，否则为对应颜色
     cex = 1)
box(lwd=3,bty="l")




adj.data <- adjust.data.by.covariates(formula=~sex+stage, data=Black_riceInfo,
                                      otu.table=otu_order, dist.method="bray")
PCs <- eigen(adj.data$adj.dist, symmetric=TRUE)
PCs.percentage = signif(PCs$values[1:2]/sum(PCs$values),3)*100

par(mfrow=c(1,1), mar=c(5, 4, 4, 8) + 0.1, pty="s", xpd=TRUE)
color = c(rep("#2E8B57",38),rep("#88ada6",15),rep("#D8BFD8",14),
          rep("#4682B4",23),rep("#808080",15),rep("#556B2F",14),
          rep("#5F9EA0",16))
plot(PCs$vectors[,1], PCs$vectors[,2],
     main="Ordination based on Bray-Curtis distance",
     xlab=paste("PC1 (", PCs.percentage[1], "%)", sep=""),
     ylab=paste("PC2 (", PCs.percentage[2], "%)", sep=""),
     col=color, pch=20)
ordiellipse(ord=PCs$vectors, groups=factor(Black_riceInfo$group, exclude=c()),
            conf=0.9, col=c("#2E8B57","#a1afc9","#88ada6","#4682B4","#D8BFD8","#808080","#800000","#556B2F","#ffc773","#5F9EA0"))
legend(x="topright", inset=c(-0.5, 0.4), bty="n",
       legend=c("WT_14w","C_14w","B_14w","W_14w","P_14w",
                "WT_22w","C_22w","B_22w","W_22w","P_22w"), pch=21,
       col=c("#2E8B57","#a1afc9","#88ada6","#4682B4","#D8BFD8","#808080","#800000","#556B2F","#ffc773","#5F9EA0"), lty=0)



################################ Co-currence #########################
#install.packages("igraph")
library(igraph)
#BiocManager::install("impute")
library(impute)
#BiocManager::install("preprocessCore")
#install.packages("WGCNA") 
library(WGCNA)


# 设置随机数种子，确保数据可重复
set.seed(123)
# 生成A、B两个样本各三个重复，共1000个ASV的丰度表
#otu <- data.frame(replicate(6, sample.int(10, 1000, replace = T))) 
#rownames(otu) <- paste0('ASV_', 1:nrow(otu)) # 行命名
#colnames(otu) <- c('A1', 'A2', 'A3', 'B1', 'B2', 'B3') # 列命名
#dim(otu) 
#View(otu)
CorrDF <- function(cormat, pmat) {
  ut <- upper.tri(cormat) # 由于相关性矩阵是对称的，取上三角矩阵计算即可
  data.frame(
    from = rownames(cormat)[col(cormat)[ut]],
    to = rownames(cormat)[row(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

connect <- function( graph , vertice, graphSize ) {
  tossResult <- toss()
  if ( tossResult == 1 ){
    nodeToConnect <- sample(1:graphSize, 1 , replace=T)
    print(nodeToConnect)
    graph <- graph + edge(vertice, nodeToConnect)
  }
  graph
}
?connect

meta <- 
  Black_riceInfo[Black_riceInfo$group%in%c("C_14w","B_14w","C_22w","B_22w"),]
w <- otu_order[,rownames(meta)]
otu_presence = which(rowSums(w)>=10)
w = w[otu_presence,]


meta22 <- Black_riceInfo[grep("22w",rownames(Black_riceInfo)),]%>%
  .[.$group%in%c("C_22w","B_22w"),]
w22 <- otu_order[,rownames(meta22)]
otu_presence = which(rowSums(w22)>=10)
w22 = w22[otu_presence,]
w22 <- w22[rownames(w22)%in%unique(diff_merge_lefse$Biomarkernames),]

occor <- corAndPvalue(t(w22),  method='spearman') # 计算ASV/OTU之间的spearman相关性
cor_df <- CorrDF(occor$cor , occor$p) # 整理ASV之间的连接关系
#View(cor_df)
#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.2),] # 保留spearman相关性绝对值>0.6的边
#cor_df <- cor_df[which(cor_df$p < 0.05),] # 保留p-value < 0.001的边
#View(cor_df)
set.seed(1)
igraph <- graph_from_data_frame(cor_df, direct=F) 
length(V(igraph)) # 查看节点数
length(E(igraph)) # 查看边数
V(igraph)
V(igraph)$size <- degree(igraph)
V(igraph)$size <- degree(igraph)*0.3# 节点的大小与节点的度成正比，进行0.8倍放缩
# 生成一些随机颜色
cols <- unique(c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                 "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                 "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                 "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                 "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                 "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f"))
V(igraph)$color <- sample(cols, length(V(igraph)), replace = T) # 从随机颜色中采样给节点上色
E(igraph)$color[E(igraph)$cor >= 0] <- "gray" # 正相关则边为深灰色
E(igraph)$color[E(igraph)$cor <= 0] <- "gray"

E(igraph)$color[E(igraph)$cor <= -0.25] <- "#3F85FF" 
E(igraph)$color[E(igraph)$cor >= 0.25] <- "#E23D30"

#E(igraph)$color[E(igraph)$cor == -0.234354982] <- "#3F85FF" 
#E(igraph)$color[E(igraph)$cor == -0.320477503] <- "#3F85FF" 
#E(igraph)$color[E(igraph)$cor == 0.421364144] <- "#E23D30"
#E(igraph)$color[E(igraph)$cor == 0.290680946] <- "#E23D30"

E(igraph)$width <- abs(E(igraph)$cor)*0.5 # 边的粗细与相关系数成正比，进行0.5倍放缩

coords <- layout_with_fr(igraph, niter=9999,grid="nogrid") # 生成布局
#set.seed(100)
plot(igraph, layout=coords)
plot(igraph, layout=coords, vertex.label =NA, vertex.frame.color=NA) # 画图

?plot


dis_bacteria<- vegan::vegdist(t(otu_order[,]), method = 'bray')
dis_archaea<- vegan::vegdist(t(archaea[,1:24]), method = 'bray')
dis_fungi<- vegan::vegdist(t(fungi[,1:24]), method = 'bray')
index<-rbind(dis_bacteria,dis_archaea)
index<-rbind(index,dis_fungi)
rownames(index)<-c("Bacteria","Archaea","Fungi")
data<-cbind(rownames(index),index)
#转化为长数据
library(reshape2)
io=melt(data,id.vars = c("V1"),value.name ="Normalized_Abundance",variable.name = "Properties")####把宽数据变为长数据
colnames(io)<-c("Gene","Treatment","Abundance")
io<-io[4:831,]
#转化为数值型
io$Abundance<-as.numeric(as.character(io$Abundance))
#绘图
ggplot(io,mapping = aes(x=Gene,y=Abundance,color=Gene))+
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.3)+
  geom_boxplot(aes(),notch = F)+
  geom_jitter(binaxis="y",position = position_jitter(0.2),satckdir="center",dotsize=0.6,size=1,alpha=0.1)+
  labs(x="", y="Beta diversity", title = "")+ theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_color_manual(values=c("#9E6C61" ,"#596E90","#599071"))+#颜色
  scale_fill_manual(values=c("#9E6C61" ,"#596E90","#599071"))+#颜色
  theme(axis.text=element_text(colour='black',size=9))


library(lme4)
install.packages("mediation")
library(mediation)
model.m=lmer(meta~microbe+age+sex+(1|time),data)
model.y=lmer(pheno~microbe+meta+age+sex+(1|time),data)
summary=summary(mediate(model.m ,model.y,treat = "microbe", mediator = "meta",boot = F,sims = 1000))

library(vegan)
library(aPCoA)
#分组和协变量数据
data <- data.frame(stage = meta$stage,sex = meta$sex,group = meta$group)
rownames(data) <- rownames(as.matrix(bray))

#执行 aPCoA 校正混杂协变量以显示 PCoA 图中主要协变量的影响，详情 ?aPCoA
#本示例中，“bray~block”意为在对群落数据进行 PCoA 时消除由位置（调查区域）带来的影响，“maincov=treatment”在作图时按试验类型分组着色
opar <- par(mfrow = c(1, 2),
            mar = c(3.1, 3.1, 3.1, 5.1),
            mgp = c(2, 0.5, 0),
            oma = c(0, 0, 0, 4))
result <- aPCoA(bray~sex+stage, data, maincov = group)
par(opar)

#提取 aPCoA 中的样本得分
aPCoA.score <- data.frame(result$plotMatrix)


dune.pairwise.adonis <- pairwise.adonis(x=aPCoA.score, 
                                        factors=meta$group, 
                                        sim.function = "vegdist",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)
dune.pairwise.adonis


########################################## Diff lefse #####################################################
diff_merge_lefse <- read.delim("F:/A_MetaboAnalysis/Black_rice/diff_merge_lefse.txt")
diff_merge_lefse$panel <- paste0(diff_merge_lefse$group,"_",diff_merge_lefse$stage)
diff_merge_lefse$group <- factor(diff_merge_lefse$group,levels=c("CvsWT","CvsB","CvsW","CvsP","PvsW","BvsW"))
diff_merge_lefse <- diff_merge_lefse[diff_merge_lefse$group%in%c("CvsWT","CvsB"),]
#View(diff_merge_lefse)
data <- otu_order[unique(diff_merge_lefse$Biomarkernames),]
#View(data)
phData <- apply(data,1,function(x){tapply(x,Black_riceInfo$group,mean)}) %>% t() %>% as.data.frame()
#View(phData)
colnames(phData)
unique(rownames(phData))
library(ComplexHeatmap)
library(circlize)
top_anno = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00")),
                       labels =c("WT_14w","C_14w","B_14w",
                                 "WT_22w","C_22w","B_22w"),
                       labels_gp = gpar(col = "white", fontsize = 14)))
select <- phData[rownames(phData)%in%c("Bifidobacterium_pseudolongum",
                                       "Bacteroides_uniformis",
                                       "Lactobacillus_johnsonii",
                                       "Escherichia_coli",
                                       "Lactobacillus_reuteri","Lactobacillus_lactis"),]


dim(select)
gene_pos <- which(rownames(phData)  %in% rownames(select))

#gene_pos <- which(rownames(phData) %in% unique(rownames(select)))
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,
                                                 labels = rownames(select)))
col_fun = colorRamp2(c(-1, 0, 1), c("cornflowerblue","white","red"))
n=t(scale(t(phData)))
?Heatmap
p1 <- Heatmap(na.omit(n),col=col_fun,border = "black",
              show_row_names = FALSE,show_column_names = F,cluster_columns = FALSE,cluster_rows = T,
              column_split = factor(colnames(phData),levels=c("WT_14w","C_14w","B_14w",
                                                              "WT_22w","C_22w","B_22w")),
              #              row_split = diff_merge_lefse$group,
              top_annotation = top_anno,
              #              left_annotation = left_anno,
              right_annotation = row_anno,row_title = NULL,column_title = NULL,
              heatmap_legend_param = list(
                title = "scale",
                title_position = "leftcenter-rot"))

p1











library(vegan)
library(ape)
library(ggplot2)
library(grid)
library(dplyr)
library(multcomp)
library(patchwork)
pcoa<- pcoa(unweight_unifrac, correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
group <- Black_riceInfo$group
plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,group)
colnames(plotdata) <-c("sample","PC1","PC2","Group")

pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)
plotdata$Group <- factor(plotdata$Group)

length=length(unique(as.character(plotdata$Group)))
times1=length%%8
times1
res1=length%%8
times2=length%%5
res2=length%%5
col1=rep(1:8,times1)
col=c(col1,1:res1)
pich1=rep(c(21:24),times2)
pich=c(pich1,15:(15+res2))
cbbPalette <- c("#B12C39","#1149A0","#E63428","#035DB3","#F65D34","#036685",
                "#F89C91","#0397C5","#F8C9B5","#87B0CC")
yf$Group
yf <- plotdata
yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1
fit1 <- aov(PC1~Group,data = plotdata)

tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
res1 <- cld(tuk1,alpah=0.05)
fit2 <- aov(PC2~Group,data = plotdata)
tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
res2 <- cld(tuk2,alpah=0.05)
test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,Group = unique(plotdata$Group))
test$Group <- factor(test$Group)

p1 <- ggplot(plotdata,aes(Group,PC1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
            size = 7,color = "black",fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=20,face = "bold"),
        axis.text.x=element_blank(),
        legend.position = "none")
p1
p3 <- ggplot(plotdata,aes(Group,PC2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
            size = 7,color = "black",fontface = "bold") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=20,angle = 45,
                                 vjust = 1,hjust = 1,face = "bold"),
        axis.text.y=element_blank(),
        legend.position = "none")

p2<-ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(fill=Group),size=8,pch = 21)+
  scale_fill_manual(values=cbbPalette,name = "Group")+
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
  ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
  theme(text=element_text(size=24))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=24),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=24),
        axis.title.y=element_text(colour='black', size=24),
        axis.text=element_text(colour='black',size=24),
        legend.title=element_text(size = 20),
        legend.text=element_text(size=16),
        legend.key=element_blank(),legend.position = c(0.95,0.25),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(0.8,"cm")) +
  guides(fill = guide_legend(ncol = 1))

p2
#?adonis
plotdata$Group%>%table()
otu.adonis=adonis(unweight_unifrac~group,data =Black_riceInfo,distance = "bray")
p4 <- ggplot(plotdata, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",
                                               otu.adonis$aov.tab$Df[1],  "\nR2 = ",
                                               round(otu.adonis$aov.tab$R2[1],4),  
                                               "\np-value = ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),size = 7) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
p5






#weight
res = pcoa(weight_unifrac,correction = "none",rn = NULL)
biplot.pcoa(res)

pcoachoose = c("Axis.1","Axis.2")
ggdata = res$vectors[,pcoachoose] %>% cbind(.,Black_riceInfo)

p = ggplot(data = ggdata,aes(x=get(pcoachoose[1]),y=get(pcoachoose[2]),colour= group))+geom_point(size= 4)+
  stat_ellipse(type = "t", linetype = 1)
mi=c("#B12C39","#1149A0","#E63428","#035DB3","#F65D34","#036685",
     "#F89C91","#0397C5","#F8C9B5","#87B0CC")
p=p+theme_classic()+scale_colour_manual(values = mi)+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text = element_text(size = 14),axis.text.x = element_text(angle = 0, hjust =0.5),
        legend.text = element_text(size=12),legend.position = "top")+xlab("PCoA1")+ylab("PCoA2")
p

#NMDS
library("amplicon", lib.loc="~/R/win-library/4.0")
library(phyloseq)
?BetaDiv
result=BetaDiv(otu=mat2, map=Black_riceInfo, group="group", 
               dist="bray", method="NMDS", Micromet="adonis")
result
#提取排序散点图(结果列表中的1)
p=result[[1]]
p
#ggsave(paste0("p3.NMDS.bray.jpg"), p, width=89, height=56, units="mm")
#ggsave(paste0("p3.NMDS.bray.pdf"), p, width=89, height=56, units="mm")

# 提取出图坐标
plotdata=result[[2]]
plotdata[1:3,1:3]

# 提取带标签排序散点图
p=result[[3]]
p$data
ggsave(paste0("p4.NMDS.bray.label.jpg"), p, width=89, height=56, units="mm")
ggsave(paste0("p4.NMDS.bray.label.pdf"), p, width=89, height=56, units="mm")

# 提取两两比较差异检测结果
pair=result[[4]]
pair
# 提取全部组整体差异检测结果
Mtest=result[[5]]
Mtest
ggplot(data=plotdata,aes(x=x,y=y))+
  geom_point(size= 4 ,aes(fill= group))+scale_colour_manual(values = cbbPalette)+
  labs(x ="",y = "")

## 稀释曲线
BiocManager::install("microbiome")
BiocManager::install("amplicon")
library(microbiome)
library(amplicon)
p = alpha_rare_curve(alpha_rare, metadata, groupID = "Group")

library(ggplot2)
library(ggdendro)
library(phyloseq)
library(tidyverse)
library(ggfun)
library(ggtree)
library(ggstance)
library(amplicon)
taxonomy <- metaphlan_anno

result <- tax_stack_clust(otu=otu,tax=taxonomy,map=Black_riceInfo,j = "Species",
                          Top=20,tran=T,dist="bray",
                          hcluter_method="ward.D",
                          cuttree=3,Group="group")
result[[3]]
?tax_stack_clust

#
tax_stack_clust <- function(
  otu=NULL,
  map=NULL,
  tax=NULL,
  dist="bray",
  Group="Group",
  j="Phylum", # 使用门水平绘制丰度图表
  rep=6 ,# 重复数量是6个
  Top=10, # 提取丰度前十的物种注释
  tran=TRUE, # 转化为相对丰度值
  hcluter_method="complete",
  cuttree=3){
  
  # 加载默认参数调试函数
  # otu=otutab
  # map=metadata
  # tax=taxonomy
  # dist="bray"
  # Group="Group"
  # j="Phylum" # 使用门水平绘制丰度图表
  # rep=6 # 重复数量是6个
  # Top=10 # 提取丰度前十的物种注释
  # tran=TRUE # 转化为相对丰度值
  # hcluter_method="complete"
  # cuttree=3
  
  # 加载包
  p_list = c("ggplot2", "BiocManager", "tidyverse", "ggdendro", "ggstance")
  for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}
  p_list = c("ggtree", "phyloseq")
  for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}
  
  # phyloseq导出特征表函数
  vegan_otu = function(physeq){
    OTU= otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU= t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  
  # phyloseq导出物种注释函数
  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)
    return(as(tax,"matrix"))
  }
  
  # 数据交叉筛选
  # Extract only those ID in common between the two tables
  idx=rownames(otu) %in% rownames(tax)
  otu=otu[idx,]
  tax=tax[rownames(otu),]
  
  # 分组列重命名为Group
  map <- map[Group]
  colnames(map)="Group"
  map$ID=row.names(map)
  
  # 数据导入phyloseq
  ps=phyloseq(sample_data(map),otu_table(as.matrix(otu), taxa_are_rows=TRUE), tax_table(as.matrix(tax)))
  # phyloseq(ps)对象标准化
  ps1_rela=phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
  
  # 导出OTU表
  otu=as.data.frame(t(vegan_otu(ps1_rela)))
  
  #计算距离矩阵
  unif=phyloseq::distance(ps1_rela , method=dist)
  # 聚类树，method默认为complete
  hc <- stats::hclust(unif, method=hcluter_method )
  
  #  take grouping with hcluster tree
  clus <- cutree(hc, cuttree)
  # 提取树中分组的标签和分组编号
  d=data.frame(label=names(clus),
               member=factor(clus))
  # Extract mapping file
  map=as.data.frame(sample_data(ps))
  # 合并树信息到样本元数据
  dd=merge(d,map,by="row.names",all=F)
  row.names(dd)=dd$Row.names
  dd$Row.names=NULL
  
  # ggtree绘图 #----
  p =ggtree(hc) %<+% dd +
    geom_tippoint(size=5, shape=21, aes(fill= Group, x=x)) +
    geom_tiplab(aes(color=Group,x=x * 1.2), hjust=1)
  # 按照分类学门(j)合并
  psdata=ps1_rela %>% tax_glom(taxrank=j)
  
  # 转化丰度值
  if (tran == TRUE) {
    psdata=psdata%>% transform_sample_counts(function(x) {x/sum(x)} )
  }
  
  #--提取otu和物种注释表格
  otu=otu_table(psdata)
  tax=tax_table(psdata)
  
  #--按照指定的Top数量进行筛选与合并
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing=TRUE)[1:Top])) {
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "Other"
    }
  }
  tax_table(psdata)= tax
  
  ##转化为表格
  Taxonomies <- psdata %>% psmelt()
  # head(Taxonomies)
  Taxonomies$Abundance=Taxonomies$Abundance * 100
  
  Taxonomies$OTU=NULL
  colnames(Taxonomies)[1]="id"
  # head(Taxonomies)
  
  p <- p + ggnewscale::new_scale_fill()
  p1 <- facet_plot(p, panel='Stacked Barplot', data=Taxonomies, geom=geom_barh,mapping=aes(x=Abundance, fill=Species),color="black",stat='identity' )
  
  grotax <- Taxonomies %>%
    group_by(Group,Species) %>%
    summarise(Abundance=mean(Abundance))
  
  #--绘制分组的聚类结果
  ps1_rela=phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
  #按照分组合并OTU表格
  hc=phyloseq::merge_samples(ps1_rela, "Group",fun=mean) %>%
    phyloseq::distance(method=dist) %>%
    stats::hclust( method=hcluter_method )
  
  #  take grouping with hcluster tree
  clus <- cutree(hc, cuttree)
  # 提取树中分组的标签和分组编号
  d=data.frame(label=names(clus), member=factor(clus))
  # Extract mapping file
  map=as.data.frame(sample_data(phyloseq::merge_samples(ps1_rela, "Group",fun=mean)))
  # 合并树信息到样本元数据
  dd=merge(d,map,by="row.names",all=F)
  row.names(dd)=dd$Row.names
  dd$Row.names=NULL
  
  # ggtree绘图 #----
  p3 =ggtree(hc) %<+% dd +
    geom_tippoint(size=5, shape=21, aes(fill= member, x=x)) +
    geom_tiplab(aes(color=member,x=x * 1.2), hjust=1)
  p3 <- p3 + ggnewscale::new_scale_fill()
  p4 <- facet_plot(p3, panel='Stacked Barplot', data=grotax, geom=geom_barh,mapping=aes(x=Abundance, fill=Species),color="black",stat='identity' )
  
  return(list(p,p1,p3,p4,Taxonomies))
}









# ABT analysis(Aggregated boosted tree（ABT）评估变量的相对重要性)
install.packages(c("dismo","gbm"))
library(dismo)
library(gbm)

#处理之前要保证两个数据文件均是行为样本。
#先计算OTU丰度表的Bray-Curtis距离。
library(vegan)
library(reshape2)
bray.otu <- vegdist(otu,method = "bray")
bray.otu <- as.matrix(bray.otu)
bray.otu <- melt(bray.otu)
colnames(z) <- colnames(env)
ABT.data <- cbind(z,bray.otu$value)
ABT.data <- as.data.frame(ABT.data)
ABT.fit <- gbm.step(data = ABT.data,gbm.x = 1:10,gbm.y = 11,family = "laplace",
                    tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)
#gbm.x的数值为1至环境因子的数目，gbm.y的数值为环境因子数目加1

#将变量重要性的结果信息提取后输出
softcorals_var_influence <- summary(ABT.fit)
write.csv(softcorals_var_influence, 'softcorals.var_influence.csv', row.names = FALSE, quote = FALSE)

#然后打开一个已经安装 ggplot2 的较新的 R 版本
#加载 ggplot2，读取数据后重新绘制变量重要性的柱形图
library(ggplot2)

softcorals_var_influence <- read.csv('softcorals.var_influence.csv', stringsAsFactors = FALSE)
softcorals_var_influence <- softcorals_var_influence[order(softcorals_var_influence$rel.influence), ]
softcorals_var_influence$var <- factor(softcorals_var_influence$var, levels = softcorals_var_influence$var)

#颜色代表不同的变量
p <- ggplot(softcorals_var_influence, aes(var, rel.influence)) +
  coord_flip() +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = '', y = 'Relative influence(%)', title = '')

p + geom_col(aes(fill = var), width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c('#8DD3C7', '#FFFFB3', '#BEBADA',
                               '#FB8072', '#80B1D3', '#FDB462', '#B3DE69'))

#或者颜色以渐变色表示变量的重要性
p + geom_col(aes(fill = rel.influence), width = 0.7, show.legend = FALSE) +
  scale_fill_gradientn(colours = colorRampPalette(c('#86CEF7', '#0000F6'))(10))

# 普氏分析结果可视化--物种与环境、物种与物种、物种与功能关联分析
# https://mp.weixin.qq.com/s/CzlW_I8fcHrHdpJ_dqVV-A
# 在微生物群落研究的过程中，我们经常需要评估微生物群落结构与环境因子整体之间是否具有显著的相关性，此时，通常使用的方式是Mantel test和普氏分析。
# 当然除了分析群落结构与环境因子的相关性之外，这两个分析还可以用于分析同一样品不同类型微生物群落之间的相关性，比如同一样品的稀有和丰富物种或者同一样品细菌和真菌群落结构的相关性。
# 不同微生物群落之间，该分析更多的还是用于分析微生物群落组成结构与其它功能基因组之间的关系，比如细菌组成结构与抗生素抗性组的相关性。
# 最后还有一种不太常用的用法，就是分析配对的两种不同类型样品微生物群落的相关性，比如河流或海洋同一位置水体和沉积物细菌群落组成结构的相关性，从而分析这两种关联样品类型中微生物群落的转移。
library(vegan)
data <- data[which(rowSums(data) > 0),]
data <- t(data)
s.dist <- vegdist(data,method = "bray")
r.dist <- vegdist(data)
mantel(s.dist,r.dist)
mantel(s.dist,r.dist,method = "spearman")
mds.s <- monoMDS(s.dist)
mds.r <- monoMDS(r.dist)
pro.s.r <- procrustes(mds.s,mds.r)
protest(mds.s,mds.r)
library(ggplot2)
Y <- cbind(data.frame(pro.s.r$Yrot), data.frame(pro.s.r$X))
X <- data.frame(pro.s.r$rotation)
Y$ID <- rownames(Y)
p <- ggplot(Y) +
  geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#B2182B", size = 1) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#56B4E9", size = 1) +
  geom_point(aes(X1, X2), fill = "#B2182B", size = 4, shape = 21) +
  geom_point(aes(MDS1, MDS2), fill = "#56B4E9", size = 4, shape = 21) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title="Correlation between community and environment") + 
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
  annotate('text', label = 'Procrustes analysis:\n    M2 = 0.8358, p-value = 0.035\nMantel test:\n    r = 0.1703, p-value = 0.04',
           x = -1.5, y = 1.2, size = 4,hjust = 0) +
  theme(plot.title = element_text(size=16,colour = "black",hjust = 0,face = "bold"))

## DMM分析
#Dirichlet Multinomial Mixtures
#Community typing with Dirichlet Multinomial Mixtures
#Dirichlet Multinomial Mixtures (DMM) 是一种用于对微生物群落分析数据进行群落分型（或聚类）的概率方法。 这是一个无限的混合模型，这意味着该方法可以推断出最佳数量的群落类型。 请注意，群落类型的数量可能会随数据大小而增长。
BiocManager::install("microbiome")
BiocManager::install("DirichletMultinomial")
library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
# Load example data
#data(dietswap)
#pseq <- dietswap
#pseq
# To speed up, only consider the core taxa
# that are prevalent at 0.1% relative abundance in 50% of the samples
# (note that this is not strictly correct as information is
# being discarded; one alternative would be to aggregate rare taxa)
#pseq.comp <- microbiome::transform(pseq, "compositional")
#head(pseq.comp)
#taxa <- core_members(pseq.comp, detection = 0.1/100, prevalence = 50/100)
#pseq <- prune_taxa(taxa, pseq)

# Pick the OTU count matrix
# and convert it into samples x taxa format
#dat <- abundances(pseq)
colnames(Black_rice_species)
count <- as.matrix(t(mat2[,1:10])) # WT_M
count <- as.matrix(t(mat2[,11:18])) # WT_F
count <- as.matrix(t(mat2[,19:25])) #W_M
count <- as.matrix(t(mat2[,26:33])) #W_F
count <- as.matrix(t(mat2[,34:40])) #P_M
count <- as.matrix(t(mat2[,41:47])) #P_F
count <- as.matrix(t(mat2[,48:55])) #C_M
count <- as.matrix(t(mat2[,56:63])) #C_F
count <- as.matrix(t(mat2[,64:70])) #B_M
count <- as.matrix(t(mat2[,71:77])) #B_F
WT <- as.matrix(t(Black_rice_species[,71:77]))

# 拟合 DMM 模型.，让我们将群落类型的最大允许数量设置为3，以加速示例。
fit <- lapply(1:3, dmn, count = count, verbose=TRUE)


lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
plot(lplc, type="b")
fit[[which.min(lplc)]]
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
#plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
#lines(aic, type="b", lty = 2)
#lines(bic, type="b", lty = 3)
#选择最佳模型

best <- fit[[which.min(unlist(lplc))]]

#参数pi及theta
mixturewt(best)
##          pi     theta
## 1 0.3738027 159.10473
## 2 0.3188891  81.91265
## 3 0.3073082  64.24696
#元素（otu）分配给不同cluster
ass <- apply(mixture(best), 1, which.max)
#每个otu对每个组成群落的贡献

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity",fill="#03989E") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))+theme_classic()+
    theme(axis.title =element_text(size = 20),
          axis.text =element_text(size = 20, color = 'black'))+
    theme(axis.text = element_text(size = 20),axis.text.x = element_text(angle = 0, hjust =0.5),
          legend.text = element_text(size=12),legend.position = "right")
  print(p)
}

## 相似性分析（ANOSIM）用来检验组间的差异是否显著大于组内差异
### ANOSIM statistic R为Anosim检验的统计量，他的分布衡量的就是零模型的分布，
### Upper quantiles of permutations就是通过999次置换获得的统计量的分位数。
### 具体说来，Anosim分析的原理是先计算样品两两之间的距离，将样品两两之间的距离按照从大到小进行排序并计算排名（秩，r），
### 并根据距离的归类（属于组间距离还是组内距离）来计算组间距离秩的均值rb与组内距离秩的均值rw之差作为统计量：
### 假如R>0，说明组内距离小于组间距离，也即分组是有效的，这与方差分析中比较组内方差与组间方差来判断的原理是类似的。
### 由上面分析结果可以看到R=0.4613，大于零模型99%分位数0.290，因此p值为0.001，结果是显著的。
# 使用 OTU 丰度表，则需要在计算时指定所依据的距离类型，这里依然使用 Bray-Curtis 距离
otu_t <- t(Black_rice_species[,1:77])%>%as.data.frame()
anosim_result_otu <- anosim(otu_t, Black_riceInfo$group, permutations = 999, distance = 'bray')        

# 首先根据丰度表计算样本距离，在将所的距离数据作为输入
bray_dis <- vegdist(otu_t, method = 'bray')
melt(bray_dis)
pheatmap(bray_dis)



anosim <- anosim(bray_dis, Black_riceInfo$group, permutations = 999)
summary(anosim)
# 可将 Dissimilarity ranks 拿出来做一个箱形图可视化
mycol=c(52,619,453,71,134,448,548,655,574,36,544,89,120,131,596,147,576)
mycol=colors()[mycol]
par(mar=c(5,5,5,5))
result=paste("R=",anosim$statistic,"p=", anosim$signif)
boxplot(anosim$dis.rank~anosim$class.vec, pch="+", col=mycol, range=1, boxwex=0.5, notch=F, ylab="Bray-Curtis Rank", main="Bray-Curtis Anosim", sub=result)

adonis <- adonis(bray_dis~group, Black_riceInfo, permutations = 999)
count1 <- as.matrix(t(cbind(mat2[,1:18],mat2[,64:77]))) # WT
# as.matrix(t(mat2[,64:77])) # Black
# as.matrix(t(mat2[,34:47])) #P
count2 <- as.matrix(t(cbind(mat2[,1:18],mat2[,34:47]))) # WT
count3 <- as.matrix(t(cbind(mat2[,64:77],mat2[,34:47]))) # WT

metadata1 <- as.data.frame(strsplit2(rownames(count1),"_"))
colnames(metadata1) <- c("group","sex","replicate")
rownames(metadata1) <- rownames(count1)
bray1<- vegdist(count1, method = 'bray')
adonis1 <- adonis(bray1~group,metadata1, permutations = 999)
adonis1

metadata2 <- as.data.frame(strsplit2(rownames(count2),"_"))
colnames(metadata2) <- c("group","sex","replicate")
rownames(metadata2) <- rownames(count1)
bray2<- vegdist(count2, method = 'bray')
adonis2 <- adonis(bray2~group,metadata2, permutations = 999)
adonis2

metadata3 <- as.data.frame(strsplit2(rownames(count3),"_"))
colnames(metadata3) <- c("group","sex","replicate")
rownames(metadata3) <- rownames(count1)
bray3<- vegdist(count3, method = 'bray')
adonis3 <- adonis(bray3~group,metadata3, permutations = 999)
adonis3

## MRPP分析
library(vegan)
single_mrpp <- function(otu_table, group_name, metadata){
  Transp_otu <- t(otu_table)
  group <- as.vector(unlist(metadata[group_name]))
  
  mrpp_all <- mrpp(Transp_otu, group, 
                   permutations = 999, 
                   distance = 'bray')
  
  results <- as.data.frame(cbind(group_name,
                                 mrpp_all$distance,
                                 round(mrpp_all$A, 3),
                                 round(mrpp_all$delta, 3),
                                 round(mrpp_all$`E.delta`, 3),
                                 round(mrpp_all$Pvalue, 3)))
  
  names(results) <- c("Group", "Distance", "Corrected_A", "Observeed_delta", 
                      "Expected_delta", "P_value")
  
  return(results)
}

mrpp_all <- single_mrpp(otu, "group", Black_riceInfo)
mrpp_all
Black_riceInfo$sample <- rownames(Black_riceInfo) 
pairwise_mrpp1 <- function(otu_table, group_name, metadata){
  if (! "sample" %in% colnames(metadata)){
    print("The colnames of your sample is not sample")
    stop()
  }else{
    group <- unique(as.vector(unlist(metadata[group_name])))
    Transp_otu <- t(otu_table)
    results <- NULL
    for (i in 1:(length(group) - 1)){
      for (j in (i + 1):length(group)){
        group_ij <- Black_riceInfo[as.vector(unlist(Black_riceInfo["condition"])) %in% c(group[i], group[j]),]
        otu_ij <- Transp_otu[group_ij$sample,  ]
        grouping <- as.vector(unlist(group_ij[group_name]))
        
        mrpp_result_ij <- mrpp(otu_ij, grouping, 
                               permutations = 999, 
                               distance = 'jaccard') 
        mrpp_result_ij
        results
        results <- rbind(results, 
                         data.frame(paste(group[i], group[j], sep = '/'), 
                                    mrpp_result_ij$distance,
                                    round(mrpp_result_ij$A, 3), 
                                    round( mrpp_result_ij$delta, 3), 
                                    round(mrpp_result_ij$`E.delta`, 3), 
                                    round( mrpp_result_ij$Pvalue, 3)))
      }
    }
  }
  results <- data.frame(results, stringsAsFactors = FALSE)
  names(results) <- c("Group", "Distance", "Corrected_A", "Observeed_delta", 
                      "Expected_delta", "P_value")
  return(results)
}
pairwise_mrpp2 <- function(otu_table, group_name, metadata){
  if (! "sample" %in% colnames(metadata)){
    print("The colnames of your sample is not sample")
    stop()
  }else{
    group <- unique(as.vector(unlist(metadata[group_name])))
    Transp_otu <- t(otu_table)
    results <- NULL
    for (i in 1:(length(group) - 1)){
      for (j in (i + 1):length(group)){
        group_ij <- subset(metadata, as.vector(unlist(metadata[group_name])) %in% c(group[i], group[j]))
        otu_ij <- Transp_otu[group_ij$sample,  ]
        grouping <- as.vector(unlist(group_ij[group_name]))
        
        mrpp_result_ij <- mrpp(otu_ij, grouping, 
                               permutations = 999, 
                               distance = 'jaccard') 
        mrpp_result_ij
        results
        results <- rbind(results, 
                         data.frame(paste(group[i], group[j], sep = '/'), 
                                    mrpp_result_ij$distance,
                                    round(mrpp_result_ij$A, 3), 
                                    round( mrpp_result_ij$delta, 3), 
                                    round(mrpp_result_ij$`E.delta`, 3), 
                                    round( mrpp_result_ij$Pvalue, 3)))
      }
    }
  }
  results <- data.frame(results, stringsAsFactors = FALSE)
  names(results) <- c("Group", "Distance", "Corrected_A", "Observeed_delta", 
                      "Expected_delta", "P_value")
  return(results)
}
#otu_t$sample <- rownames(otu_t)
#rownames(otu_t) <- NULL
#otu_t <- otu_t %>%
#  dplyr::select(sample, everything())
#colnames(otu_t)
head(Black_riceInfo)
pairwise_results <- pairwise_mrpp1(otu, "condition", Black_riceInfo)
pairwise_results
## Significance of delta("P_value)，即显著性p值，小于0.05说明差异显著；
## Chance corrected within-group agreement A(Corrected_A)，即A值，大于0说明组间差异大于组内差异，
##小于0说明组内差异大于组间差异；observe delta值越小说明组内差异小，expect delta值越大说明组间差异大。
#Group Distance Corrected_A Observeed_delta Expected_delta P_value
#1   WT/W     bray       0.091           0.498          0.549   0.001
#2   WT/P     bray       0.111           0.467          0.526   0.001
#3   WT/C     bray       0.036           0.494          0.512   0.002
#4   WT/B     bray       0.116           0.466          0.527   0.001
#5    W/P     bray       0.032           0.486          0.502   0.009
#6    W/C     bray       0.074           0.514          0.555   0.001
#7    W/B     bray       0.044           0.485          0.507   0.001
#8    P/C     bray       0.080           0.481          0.523   0.001
#9    P/B     bray       0.090           0.448          0.493   0.001
#10   C/B     bray       0.123           0.480          0.547   0.001
pairwise_results <- pairwise_mrpp2(otu, "group", Black_riceInfo)
pairwise_results

# test
group <- unique(as.vector(unlist(Black_riceInfo["condition"])))
Transp_otu <- t(otu)
results <- NUL

i=1
j=2

for (i in 1:(length(group) - 1)){
  for (j in (i + 1):length(group)){
    group_ij <- Black_riceInfo[as.vector(unlist(Black_riceInfo["condition"])) %in% c(group[i], group[j]),]
    otu_ij <- Transp_otu[group_ij$sample,  ]
    grouping <- as.vector(unlist(group_ij["condition"]))
    
    mrpp_result_ij <- mrpp(otu_ij, grouping, 
                           permutations = 999, 
                           distance = 'bray') 
    mrpp_result_ij
    results
    results <- rbind(results, 
                     data.frame(paste(group[i], group[j], sep = '/'), 
                                mrpp_result_ij$distance,
                                round(mrpp_result_ij$A, 3), 
                                round( mrpp_result_ij$delta, 3), 
                                round(mrpp_result_ij$`E.delta`, 3), 
                                round( mrpp_result_ij$Pvalue, 3)))
    
  }}
results <- data.frame(results, stringsAsFactors = FALSE)
names(results) <- c("Group", "Distance", "Corrected_A", "Observeed_delta", 
                    "Expected_delta", "P_value")
return(results)


# lefse分析


## 
BiocManager::install("ggdendro")
BiocManager::install("ggtree")
library(ggplot2)
library(ggdendro)
library(phyloseq)
library(tidyverse)
library(ggfun)
library(ggtree)
library(ggstance)
library(amplicon)
taxonomy <- metaphlan_anno
result <- tax_stack_clust(otu=otu,tax=taxonomy,map=Black_riceInfo,
                          Top=10,tran=F,dist="bray",
                          hcluter_method="ward.D2",
                          cuttree=3,Group="group")
result[[4]]

meta <- Black_riceInfo
# general
cat('Starting confounder testing script\n')
start.time <- proc.time()[1]
set.seed(2022)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("The analysis tag needs to be provided! Exiting...\n")
}
tag <- args[1]
alpha.meta <- 1e-05

# Get Data

all_rel.red <- otu_order


# preprocess confounder variables to test later
meta$num <- 1:length(meta$condition)
meta$state <- 0
meta$state[1:52] <- "healthy"
meta$state[53:179] <- "CRC"
meta$diet <- 0
meta$diet[1:90] <- "Control"
meta$diet[91:119] <- "Black rice"
meta$diet[120:149] <- "Brown rice"
meta$diet[150:179] <- "Polish rice"



####
#  variance explained by disease status
ss.disease <- apply(all_rel.red, 1, FUN=function(x, label){
  rank.x <- rank(x)/length(x)#百分比排名
  ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)#方差
  ss.o.i <- sum(vapply(unique(label), function(l){
    sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
  }, FUN.VALUE = double(1)))/length(rank.x)#方差
  return(1-ss.o.i/ss.tot)
}, label=meta %>% pull(diet))

# calculate trimmed mean abundance
t.mean <- apply(all_rel.red, 1, mean, trim=0.1)

df.plot.all <- tibble(
  species=rownames(all_rel.red),
  disease=ss.disease,
  t.mean=t.mean)

###
# Test all possible confounder variables
df.list <- list()
for (meta.var in c('stage', 'sex', "state", "diet")){
  
  cat('###############################\n', meta.var, '\n')
  meta.c <- meta %>%
    filter(!is.na(eval(parse(text=meta.var))))
  
  cat('After filtering, the distribution of variables is:\n')
  print(table(meta.c$sex, meta.c %>% pull(meta.var)))
  
  feat.red <- all_rel.red
  
  cat('Calculating variance explained by meta-variable...\n')
  ss.var <- apply(feat.red, 1, FUN=function(x, label){
    rank.x <- rank(x)/length(x)
    ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
    ss.o.i <- sum(vapply(unique(label), function(l){
      sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
    return(1 - ss.o.i/ss.tot)
  }, label=meta.c %>% pull(meta.var))
  df.plot.all[[meta.var]] <- ss.var
  
  cat('Calculating association with the meta-variable...\n')
  if (meta.c %>% pull(meta.var) %>% unique %>% length > 2){
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      kruskal.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  } else {
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      wilcox.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  }
  meta.significance <- p.adjust(meta.significance, method='fdr')
  df.plot.all[[paste0(meta.var, '.significance')]] <- meta.significance
  cat('\n')
}

############################### 3.1 plot###############################################
# data_text<-data.frame(label=c("r=-0.013 ","r=0.18","r=1","r=0.10"),
#                       facet=c('Age','Sex','State','Study'),
#                       x=c(0.05,0.05,0.05,0.05),
#                       y=c(0.2,0.2,0.2,0.2))
f_labels <- data.frame(facet=c("diet",'sex','stage',"state"))

f_label =c("r=-0.013 ","r=0.18","r=1","r=0.10")



#使用geom_text语法映射分面的相同位置不同的标签


library(tidyverse)
# plot
library(ggplot2)
library(scales)
colnames(df.plot.all)

data=df.plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) 




p2=ggplot(data, aes(x=disease, y=meta, size=t.mean+1e-08)) +
  geom_point(shape=19) +
  
  xlab('Variance in diets') +
  ylab('Variance in different factors') +
  theme_bw() +
  facet_wrap(~type, ncol=2) +
  theme(
    axis.text = element_text(size=15, color="black"),
    #text = element_text(size=15),
    axis.text.x=element_text(hjust=0.5, size=25,face = "bold"),
    axis.title = element_text(size=28,face = "bold"),
    axis.text.y = element_text(size =25,face = "bold"),
    legend.text=element_text(size=25,face = "bold"),
    legend.title=element_text(size=25,face = "bold"),
    strip.text = element_text(size=30,face="bold",color="black"),    #控制分面标题的文字
    strip.background = element_blank(),###分页标题背景
    panel.grid.minor = element_blank()##网格线
  )+
  
  scale_x_continuous(breaks = seq(from=0, to=0.6, by=0.1)) +
  scale_y_continuous(breaks=seq(from=0, to=0.6, by=0.1))+
  scale_size_area(name='Trimmed mean\nabundance',
                  breaks=c(1e-05, 1e-03, 1e-02)) +
  guides( size = "legend", colour='legend')
#geom_text(x = 0.05, y =0.2, aes(label = f_label), data = f_labels)
#geom_text(aes(x=x,y=y, label=label,colour=NULL),data=data_text)+
p2


ggsave(p2,filename ='/home/panxl/CRC_Metagenomic_Metab/fig/Confounding_factors_plot.pdf',
       width = 9, height = 7, useDingbats=FALSE)




df.plot.study1=as.data.frame(df.plot.study)

cor(df.plot.all$disease ,df.plot.all$state)#1
cor(df.plot.all$disease ,df.plot.all$stage)#0.1270992
cor(df.plot.all$disease ,df.plot.all$sex)#-0.1790837
cor(df.plot.all$disease ,df.plot.all$diet)# 0.5315652

trpB <- read.delim("F:/A_MetaboAnalysis/Black_rice/trpB.txt")
a <- strsplit2(trpB$Group,"_")
colnames(a) <- c("Condition","Stage")
trpB <- cbind(trpB,a)
trpB$Group <- factor(trpB$Group,levels = c("WT_14w","C_14w","B_14w",
                                           "WT_22w","C_22w","B_22w"))

ggplot(trpB,aes(Group,log2(TrpB+1)))+
  stat_summary(mapping=aes(fill = Group),fun=mean,geom = "bar",
               fun.args = list(mult=1),width=0.5)+
  stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2,lwd=1)+
  theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+
  stat_summary(fun=mean, geom="point", size=2)+ylab("TrpB relative abundance")+xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')



ggplot(trpB, aes(x=Group, y=TrpB, fill=Group)) +
  geom_jitter() +
  scale_colour_ordinal()+
  # geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=0.6,position=position_dodge(0.9))+
  theme(legend.position="none") +
  scale_fill_manual(values=c("#0ECBC5","#FFDD93","#E56A54","#0293A1","#F5C31E","#F73E00"))+						
  ylab("TrpB relative abundance") +theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold",angle = 45,hjust = 1),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+xlab("")+
  geom_signif(comparisons = list(c("WT_14w", "C_14w"),c("WT_22w","C_22w"),
                                 c("C_14w","B_14w"),c("C_22w","B_22w"),
                                 c("WT_14w", "B_14w"),c("WT_22w","B_22w")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)

