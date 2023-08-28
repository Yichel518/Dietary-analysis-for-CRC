setwd("F:/A_MetaboAnalysis/Black_rice/FMT/")
BR_newcount <- read.delim("F:/A_MetaboAnalysis/Black_rice/FMT/new_count.txt", row.names=1)
colnames(BR_newcount)
BR_newcount$B4 <- NULL
BR_newcount$C4 <- NULL
#BiocManager::install("rtracklayer")
library(rtracklayer)
library(tidyverse)
library(limma)
library(DESeq2)
library(reshape2)
library(dplyr)
library(org.Mm.eg.db)
gtf = rtracklayer::import("F:/A_MetaboAnalysis/Black_rice/weeks_22/Mus_musculus.GRCm39.104.gtf")
gtf = as.data.frame(gtf)
tra = gtf[gtf$type=="transcript",
          c("start","end","gene_id")]
glt = mutate(tra, length = end - start) %>%
  .[order(.$length),] %>% 
  filter(!duplicated(gene_id)) 
final_gene <- as.data.frame(glt[,3:4])
rownames(final_gene) <- final_gene$gene_id
final_gene <- final_gene[rownames(BR_newcount),]
kb <- final_gene$length / 1000
rpk <- BR_newcount/ kb
TPM <- t(t(rpk)/colSums(rpk) * 1000000)%>%as.data.frame()
final_gene <- final_gene[rownames(BR_newcount),]
kb <- final_gene$length / 1000
rpk <- BR_newcount/ kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)%>%as.data.frame()
pheatmap::pheatmap(cor(tpm,method = "spearman"),display_numbers = T,cluster_cols = F,cluster_rows = F)

library(PCAtools)
pca <- pca(tpm, metadata =coldata)
biplot(pca, x = 'PC1', y = 'PC2') #可以看到这两个成分对样品的解释度
#screen(p) #可以看到所有成分对样品的解释
##将pca与关联样本信息表 
pca$variance
pca_rotated_plus <- cbind(pca$rotated,coldata)
library(ggsci)
pca_rotated_plus$group <- factor(pca_rotated_plus$group,levels=c("C","B"))

ggplot(pca_rotated_plus,aes(x = PC1 , y = PC2))+
  geom_point(size = 8,aes(#shape = group, 
    fill = condition,color = condition)) +
  stat_ellipse(aes(color = condition,fill=condition),linetype = 'dashed',size = 1, show.legend = FALSE) + #添加分组椭圆
  labs(x = 'PC1 (23.72%)',y = 'PC2 (14.89%)') + 
  scale_fill_manual(values = c("#D14424","#0094A5"))+
  scale_color_manual(values = c("#D14424","#0094A5"))+
  theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.position = "top",legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))












colnames(BR_newcount)
data <- BR_newcount
coldata <- substr(colnames(BR_newcount),1,1)%>%as.data.frame(.)
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data,
                             colData = coldata,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)

res <- results(dds2,contrast=list(c("condition_C_vs_B")))%>%as.data.frame()
res$ENSEMBL <- rownames(res)
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=res$ENSEMBL,				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
res <- merge(res,sym,all.x=T,by="ENSEMBL")

meta <- res[res$SYMBOL%in%c("Cxcl10","Ccl17","Ccl2","Ccl7","Cdkn1a","Gabra3","Oprk1","Bche","Cnr1"),]
meta

res$type <- ifelse(res$pvalue < 0.05,
                   ifelse(abs(res$log2FoldChange) > 0 ,
                          ifelse(res$log2FoldChange < 0 ,'down','up'),'noSig'),'noSig')


table(res$type)
p=ggplot(na.omit(res),aes(x = log2FoldChange,y = -log10(pvalue))) +
  geom_point(aes(color = type),size = 3) +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#DA1212','noSig'='grey','down'='#3E7C17'),
                     # legend labels
                     label = c('up'='up (num=1606)','noSig'='noSig (num=33435)','down'='down (num=1241)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('GF-Control vs Black')+theme_classic()
p + geom_text_repel(data = meta,aes(x = log2FoldChange,y = -log10(pvalue),label = meta$SYMBOL),
                    force=20,color="black",size=6,point.padding = 0.5,hjust = 0.5,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                    segment.color="black",segment.size=0.2,nudge_y=1)+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))



C <- subset(res, pvalue< 0.05 & log2FoldChange > 0 )
C <- as.data.frame(C)
B <- subset(res, pvalue < 0.05 & log2FoldChange < 0 )
B <- as.data.frame(B)
C$condition <- "C"
B$condition <- "B"
Diff_CB <- rbind(C,B)
View(Diff_CB)
library(clusterProfiler)


BP1<- clusterProfiler::simplify(enrichGO(gene=rownames(C),keyType = "ENSEMBL",
                                      OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                      pvalueCutoff=1,qvalueCutoff=1,readable=T))%>%as.data.frame()


BP2<- clusterProfiler::simplify(enrichGO(gene=rownames(B),keyType = "ENSEMBL",
                                      OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                      pvalueCutoff=1,qvalueCutoff=1,readable=T))%>%as.data.frame()

BP1 <- subset(BP1,pvalue<0.05)
BP2 <- subset(BP2,pvalue<0.05)
BP1<- clusterProfiler::simplify(enrichGO(gene=C$ENSEMBL,keyType = "ENSEMBL",
                                         OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                         pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T))%>%as.data.frame()


BP2<- clusterProfiler::simplify(enrichGO(gene=B$ENSEMBL,keyType = "ENSEMBL",
                                         OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                         pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T))%>%as.data.frame()


View(BP1)
View(BP2)

write.table(BP1,"BP_C.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(BP2,"BP_B.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(BP1,"BP_C_adj.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(BP2,"BP_B_adj.txt", sep = '\t', col.names = NA, quote = FALSE)

select <- read.delim("select.txt")
select$RichFactor=as.numeric(select$Count)/ as.numeric(sub("/\\d+", "",select$BgRatio))						
select$group <- factor(select$Class,levels = c("GF-Control","GF-Black"))
ggplot(select,aes(group,Description)) + geom_point(aes(color=pvalue,size=RichFactor))+						
  scale_colour_gradient(low="#FF603F",high= "#6997ED")+						
  labs(size="RichFactor",x="",y="",title="")+					
  theme(axis.ticks = element_blank())+
  scale_y_discrete(limits=rev(unique(select$Description)))+theme_classic()+	
  theme(axis.title =element_text(size = 20,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 20,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.title  = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))

library(clusterProfiler)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(C),
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")

KEGG1<- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",
                                    organism ="mmu",
                                    pvalueCutoff = 1,qvalueCutoff = 1)%>%as.data.frame()
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(B),
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")
View(sym)
KEGG2<- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",
                                    organism = 'mmu',
                                    pvalueCutoff = 1,qvalueCutoff = 1)%>%as.data.frame()

KEGG1 <- subset(KEGG1,pvalue<0.05)
KEGG2 <- subset(KEGG2,pvalue<0.05)

View(KEGG1)
View(KEGG2) 
a <- unique(c("Tigit","Cxcl10","Ccl17","Ccl2","Ccl7","Cdkn1a","Gabra3","Oprk1","Bche","Cnr1","Cxcl16"))
library(reshape2)
library(dplyr)

tpm$ENSEMBL <- rownames(tpm)
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=tpm$ENSEMBL,				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
tpm <- merge(tpm,sym,all.x=T,by="ENSEMBL")
tpm <- na.omit(tpm[!duplicated(tpm$SYMBOL),])
rownames(tpm) <- tpm$SYMBOL
tpm$ENSEMBL <- NULL
tpm$ENTREZID <- NULL
tpm$SYMBOL <- NULL
tmp <- tpm[rownames(tpm)%in%a,]
tmp$gene <- rownames(tmp)
View(tmp)
ggdata <- melt(tmp)
ggdata$group <- factor(substr(ggdata$variable,1,1),levels = c("C","B"))
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
#Data_summary <- summarySE(ggdata, measurevar="value", groupvars=c("group","gene"))
library(ggsignif)
data <- ggdata[ggdata$gene%in%c("Ccl2"),]
a <- unique(c("Tigit","Cxcl10","Ccl17","Ccl2","Ccl7","Cdkn1a","Gabra3","Oprk1","Bche","Cnr1","Cxcl16"))

ggplot(data,aes(x = group,y = log2(value+1),fill=group))+
  stat_boxplot(geom = "errorbar", width=0.1,lwd=1.2)+
  geom_boxplot(notch = F,lwd=1.2,fatten = 1.2,width=0.5)+theme_classic()+
  scale_fill_manual(values=c("#EBBE70","#C3B3D5"))+
  geom_signif(comparisons = list(c("C", "B")),
              map_signif_level = F,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative expression (TPM)")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  theme(legend.position = 'none')+ggtitle("Ccl2")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))
