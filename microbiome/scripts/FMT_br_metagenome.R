setwd("F:/A_MetaboAnalysis/Black_rice/FMT/")
#install.packages("ggiraph")
#install.packages("rlang")
#install.packages("vctrs")
library(limma)
library(reshape2)
library(dplyr)
# 种 species
FMT <- read.delim("F:/A_MetaboAnalysis/Black_rice/FMT/merged_abundance_table_species.txt",row.names = 1 ,header=T)

library(reshape2)
library(dplyr)
library(tidyverse)
library(limma)
colnames(FMT)
dim(FMT)
Black_riceInfo <- strsplit2(colnames(FMT),"_")[,1]%>%as.data.frame()
Black_riceInfo$V2 <- substr(Black_riceInfo$.,1,1)
colnames(Black_riceInfo) <- c("sample","group")
rownames(Black_riceInfo) <- colnames(FMT)
Black_riceInfo <- Black_riceInfo[order(Black_riceInfo$group,decreasing = F),]
FMT <- FMT[,rownames(Black_riceInfo)]
colnames(FMT)
# 通常我们只关注高丰度且显著差异的，按每个OTU的总丰度排序
FMT$sum <- rowSums(FMT)
# 按每行的和降序排列
FMT_order <- FMT[order(FMT$sum, decreasing=TRUE), ]
otu <- FMT_order
otu$sum <- rowSums(otu)
otu_order <- otu[order(otu$sum, decreasing=TRUE), ]
colnames(otu)
otu_order$sum <- NULL
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

count<- read.delim("F:/A_MetaboAnalysis/Black_rice/FMT/merged_count_species.txt", row.names=1)
colnames(count)
dim(count)
count <- count[,order(colnames(count))]
meta <- substr(colnames(count),1,1)%>%as.data.frame()
colnames(meta) <- c("group")
rownames(meta) <- colnames(count)
otu_t <- t(count)
rownames(otu_t)

alpha <- alpha_diversity(otu_t)
alpha$ACE<-NULL
alpha$goods_Coverage <- NULL

alpha$group <- meta$group

head(alpha)
otu_diversity_long <- melt(alpha,id.vars = c("group"))
otu_diversity_long$group <- factor(otu_diversity_long$group,levels=c("C","B"))
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


Data_summary <- summarySE(otu_diversity_long, measurevar="value", groupvars=c("group"))
library(rstatix)
#data <- otu_diversity_long[otu_diversity_long$variable=="Shannon",]
#data$num <- 1:length(data$group)
#data[c(33:46,63:76),]%>%wilcox_test(value~group|sex, data=.)
#model <- lm(value~group, data = data)
#summary(model)
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


ggplot(data,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.1,lwd=1.2)+
  geom_boxplot(notch = F,lwd=1.2,fatten = 1.2,width=0.5)+theme_classic()+
  scale_fill_manual(values=c("#EBBE70","#B4CEE3"))+
# geom_signif(comparisons = list(c("C","B")),
#              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
#              tip_length = 0,vjust = 0.2)+
  ylab("Simpson diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 0,vjust = 1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')


data <- otu_diversity_long[otu_diversity_long$variable=="Shannon",]
data1 <- data %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup() 

ggplot(data1,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#EBBE70","#B4CEE3"))+
  geom_signif(comparisons = list(c("C","B")),
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

ggplot(data1,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#EBBE70","#B4CEE3"))+
  geom_signif(comparisons = list(c("C","B")),
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


ggplot(otu_diversity_long[otu_diversity_long$variable=="Chao1",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#EBBE70","#B4CEE3"))+
  geom_signif(comparisons = list(c("C","B")),
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
                                           "Lactococcus_lactis",
                                           "Bifidobacterium_animalis"),]
data$species <- rownames(data)
ggdata <- melt(data)
ggdata
ggdata$group <- substr(ggdata$variable,1,1)%>%factor(.,levels = c("C","B"))
data <- ggdata[ggdata$species=="Bacteroides_uniformis",]

data1 <- data %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()   
ggplot(na.omit(data1),aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.1,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2,width=0.5)+theme_classic()+
  scale_fill_manual(values=c("#EBBE70","#B4CEE3"))+
  geom_signif(comparisons = list(c("C","B")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 0, vjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  ylab("Relative abundance")+ggtitle("B.uniformis")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))


data <- ggdata[ggdata$species=="Bifidobacterium_animalis",]
data1 <- data %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()   
ggplot(na.omit(data1),aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.1,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2,width=0.5)+theme_classic()+
  scale_fill_manual(values=c("#EBBE70","#B4CEE3"))+
  geom_signif(comparisons = list(c("C","B")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 0, vjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  ylab("Relative abundance")+ggtitle("B.animalis")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))


####################################beita#########################################
?pcoa
library(ape)
metaphlan_anno = read.table("F:/A_MetaboAnalysis/metaphlan_anno.txt",header = T) 
metaphlan_tree = ape::read.tree("F:/A_MetaboAnalysis/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk")
metaphlan_tree$tip.label <- metaphlan_tree$tip.label %>% gsub(".+\\|s__","",.)
my_tree = ape::keep.tip(metaphlan_tree,intersect(rownames(otu_order),metaphlan_tree$tip.label))
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
meta <- substr(colnames(count),1,1)%>%as.data.frame()
head(meta)
colnames(meta) <-"group" 
rownames(meta) <- colnames(count)
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
ggdata$group <- factor(meta$group,levels = c("C","B"))
ggplot(data = ggdata,aes(x=get(pcoachoose[1]),y=get(pcoachoose[2]),fill=group,color = group))+						
  geom_point(size=8) +
    stat_ellipse(type = "t", linetype = 1, show.legend = FALSE,aes(color = group))+						
  ylab("Bray-Curtis diversity") +theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+						
  scale_color_manual(values=c("#EBBE70","#B4CEE3"))+						
  scale_fill_manual(values=c("#EBBE70","#B4CEE3"))+						
  xlab("PCoA1 (42.2%)")+ylab("PCoA2 (25.9%)")


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

