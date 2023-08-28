setwd("F:/A_MetaboAnalysis/Black_rice/weeks_22/")
WT_Apc_count <- read.delim("F:/A_MetaboAnalysis/Black_rice/weeks_22/WT_Apc_count.txt", row.names=1)
colnames(WT_Apc_count)
BR_newcount <- read.delim("F:/A_MetaboAnalysis/Black_rice/weeks_22/BR_newcount.txt", row.names=1)
#BiocManager::install("rtracklayer")
library(rtracklayer)
library(tidyverse)
gtf = rtracklayer::import("./Mus_musculus.GRCm39.104.gtf")
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
tpm <- t(t(rpk)/colSums(rpk) * 1000000)%>%as.data.frame()

pheatmap::pheatmap(cor(tpm,method = "spearman"),
                   color = colorRampPalette(c("white", "firebrick3"))(50))

coldata <- strsplit2(colnames(WT_Apc_count),"_")%>%as.data.frame(.)
rownames(coldata) <- colnames(WT_Apc_count)
colnames(coldata) <- c("condition","sex","tissue","replicate")
coldata$group <- paste0(coldata$condition,"_",coldata$sex,"_",coldata$tissue)
dds<- DESeqDataSetFromMatrix(countData =                                                                
                                 WT_Apc_count,
                                 colData = coldata,
                                 design = ~sex+condition+tissue)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("groupc2vsc1")))
vsd <- vst(dds2)
data <- assay(vsd)
res <- as.data.frame(res)
res$gene <- rownames(res)
summary(res)
#View(res)


#####################WT vs  Apc Cecum##################################
colnames(WT_Apc_count)
data <- WT_Apc_count[,c(1:3,7:9,13:15,19:21)]
coldata <- strsplit2(colnames(data),"_")%>%as.data.frame(.)
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data,
                             colData = coldata,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_WT_vs_C")))
WT_Cecum <- subset(res, padj < 0.05 & log2FoldChange > 0 )
WT_Cecum <- as.data.frame(WT_Cecum)
Apc_Cecum <- subset(res, padj < 0.05 & log2FoldChange < 0 )
Apc_Cecum <- as.data.frame(Apc_Cecum)

#####################WT vs  Apc Colon##################################
colnames(WT_Apc_count)
data <- WT_Apc_count[,c(4:6,10:12,16:18,22:24)]
colnames(data)
coldata <- strsplit2(colnames(data),"_")%>%as.data.frame(.)
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data,
                             colData = coldata,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_WT_vs_C")))%>%as.data.frame()
res$ENSEMBL <- rownames(res)
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=res$ENSEMBL,				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
tpm$ENSEMBL <- rownames(tpm)
res <- merge(res,sym,all.x=T,by="ENSEMBL")
View(res)
dd <- tpm[,grep("Colon",colnames(tpm))]
dd$ENSEMBL <- rownames(dd)
res <- merge(dd,res,all.x=T,by="ENSEMBL")
write.table(res,"res_WTvsApc.txt", sep = '\t', col.names = NA, quote = FALSE)

WT_Colon <- subset(res, pvalue < 0.05 & log2FoldChange > 0 )
WT_Colon <- as.data.frame(WT_Colon)
WT_Colon$Class <- "WT"
Apc_Colon <- subset(res, pvalue < 0.05 & log2FoldChange < 0 )
Apc_Colon <- as.data.frame(Apc_Colon)
Apc_Colon$Class <- "CRC"
WTvsCRC_Diff <- rbind(WT_Colon,Apc_Colon)
write.table(WTvsCRC_Diff,"WTvsCRC_Diff.txt", sep = '\t', col.names = NA, quote = FALSE)
View(res)

BP<- as.data.frame(simplify(enrichGO(gene=Apc_Colon$ENSEMBL,keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                     pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)))
BP<- as.data.frame(simplify(enrichGO(gene=WT_Colon$ENSEMBL,keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                     pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)))

View(BP)
meta <- res[res$SYMBOL%in%c("Slc11a2","Lox","Alas2","Icos","Ctla4",
                            "S100a9","Tigit","Lrrc32","Il17a","Lcn2",
                            "Cxcl12","Tnfsf14","Efnb2","Gabrg2","Cbln2",
                            "Nptx1","Add2","Adgrb3","Nrxn2","Hoxc11"),]
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
                     label = c('up'='up (num=625)','noSig'='noSig (num=34694)','down'='down (num=588)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('WT vs Apc')+theme_classic()
p + geom_text_repel(data = meta,aes(x = log2FoldChange,y = -log10(pvalue),label = meta$SYMBOL),
                    force=20,color="black",size=6,point.padding = 0.5,hjust = 0.5,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                    segment.color="black",segment.size=0.2,nudge_y=1)+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))








#####################WT M vs F Cecum##################################
colnames(WT_Apc_count)
data1 <- WT_Apc_count[,c(1:3,7:9)]
coldata1 <- strsplit2(colnames(data1),"_")%>%as.data.frame(.)
rownames(coldata1) <- colnames(data1)
colnames(coldata1) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                             data1,
                             colData = coldata1,
                             design = ~sex)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("sex_M_vs_F")))
WTM_Cecum <- subset(res, padj < 0.05 & log2FoldChange > 0 )
WTM_Cecum <- as.data.frame(WTM_Cecum)
WTF_Cecum <- subset(res, padj < 0.05 & log2FoldChange < 0 )
WTF_Cecum <- as.data.frame(WTF_Cecum)
dim(WTF_Cecum)
dim(WTM_Cecum)
library(clusterProfiler)
BP<- as.data.frame(simplify(enrichGO(gene=rownames(WTF_Cecum),keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                     pvalueCutoff=1,qvalueCutoff=1,readable=T)))
View(BP)

###############WT M vs F Colon################
colnames(WT_Apc_count)
data2 <- WT_Apc_count[,c(4:6,10:12)]
coldata2 <- strsplit2(colnames(data2),"_")%>%as.data.frame(.)
rownames(coldata2) <- colnames(data2)
colnames(coldata2) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data2,
                             colData = coldata2,
                             design = ~sex)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("sex_M_vs_F")))
WTMF_Colon <- subset(res, padj < 0.05 & log2FoldChange > 0 )
WTMF_Colon <- as.data.frame(WTMF_Colon)
WTFM_Colon <- subset(res, padj < 0.05 & log2FoldChange < 0 )
WTFM_Colon <- as.data.frame(WTFM_Colon)
write.table(WTMF_Colon,"WTMF_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(WTFM_Colon,"WTFM_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)

dim(WTF_Colon)
dim(WTM_Colon)
BP1<- as.data.frame(simplify(enrichGO(gene=rownames(WTF_Colon),keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                     pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)))
View(BP1)

BP2<- as.data.frame(simplify(enrichGO(gene=rownames(WTM_Colon),keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                     pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)))
View(BP2)
########################WTM vs ApcM Cecum###########################
colnames(WT_Apc_count)
data3 <- WT_Apc_count[,c(1:3,13:15)]
coldata3 <- strsplit2(colnames(data3),"_")%>%as.data.frame(.)
rownames(coldata3) <- colnames(data3)
colnames(coldata3) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data3,
                             colData = coldata3,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_WT_vs_C")))
WTM_Cecum <- subset(res, padj < 0.05 & log2FoldChange > 0 )
WTM_Cecum <- as.data.frame(WTM_Cecum)
ApcM_Cecum <- subset(res, padj < 0.05 & log2FoldChange < 0 )
ApcM_Cecum <- as.data.frame(ApcM_Cecum)
dim(WTM_Cecum)
dim(ApcM_Cecum)
BP<- as.data.frame(simplify(enrichGO(gene=rownames(ApcM_Cecum),keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                     pvalueCutoff=1,qvalueCutoff=1,readable=T)))
View(BP)
###########################WT M vs Apc M Colon###########################
colnames(WT_Apc_count)
data4 <- WT_Apc_count[,c(4:6,16:18)]
coldata4 <- strsplit2(colnames(data4),"_")%>%as.data.frame(.)
rownames(coldata4) <- colnames(data4)
colnames(coldata4) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data4,
                             colData = coldata4,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_WT_vs_C")))
WT_CM_Colon <- subset(res, padj < 0.05 & log2FoldChange > 0 )
WT_CM_Colon <- as.data.frame(WT_CM_Colon)
C_WTM_Colon <- subset(res, padj < 0.05 & log2FoldChange < 0 )
C_WTM_Colon <- as.data.frame(C_WTM_Colon)
write.table(WT_CM_Colon,"WT_CM_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(C_WTM_Colon,"C_WTM_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)

########################WT F vs Apc F Cecum###########################
colnames(WT_Apc_count)
data5 <- WT_Apc_count[,c(7:9,19:21)]
coldata5 <- strsplit2(colnames(data5),"_")%>%as.data.frame(.)
rownames(coldata5) <- colnames(data5)
colnames(coldata5) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data5,
                             colData = coldata5,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_WT_vs_C")))
WTF_Cecum <- subset(res, padj < 0.05 & log2FoldChange > 0 )
WTF_Cecum <- as.data.frame(WTF_Cecum)
ApcF_Cecum <- subset(res, padj < 0.05 & log2FoldChange < 0 )
ApcF_Cecum <- as.data.frame(ApcF_Cecum)
dim(WTF_Cecum)
dim(ApcF_Cecum)
###########################WT F vs Apc F Colon###########################
colnames(WT_Apc_count)
data6 <- WT_Apc_count[,c(10:12,22:24)]
coldata6 <- strsplit2(colnames(data6),"_")%>%as.data.frame(.)
rownames(coldata6) <- colnames(data6)
colnames(coldata6) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data6,
                             colData = coldata6,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_WT_vs_C")))
WT_CF_Colon <- subset(res, padj < 0.05 & log2FoldChange > 0 )
WT_CF_Colon <- as.data.frame(WT_CF_Colon)
C_WTF_Colon <- subset(res, padj < 0.05 & log2FoldChange < 0 )
C_WTF_Colon <- as.data.frame(C_WTF_Colon)
write.table(WT_CF_Colon,"WT_CF_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(C_WTF_Colon,"C_WTF_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)

########################### M vs Apc F Cecum###########################
colnames(WT_Apc_count)
data7 <- WT_Apc_count[,c(13:15,19:21)]
coldata7 <- strsplit2(colnames(data7),"_")%>%as.data.frame(.)
rownames(coldata7) <- colnames(data7)
colnames(coldata7) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data7,
                             colData = coldata7,
                             design = ~sex)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("sex_M_vs_F")))
ApcM_Cecum <- subset(res, padj < 0.05 & log2FoldChange > 0 )
ApcM_Cecum <- as.data.frame(ApcM_Cecum)
ApcF_Cecum <- subset(res, padj < 0.05 & log2FoldChange < 0 )
ApcF_Cecum <- as.data.frame(ApcF_Cecum)
dim(ApcF_Cecum)
dim(ApcM_Cecum)
########################### M vs Apc F Colon###########################
colnames(WT_Apc_count)
data8 <- WT_Apc_count[,c(16:18,22:24)]
coldata8 <- strsplit2(colnames(data8),"_")%>%as.data.frame(.)
rownames(coldata8) <- colnames(data8)
colnames(coldata8) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data8,
                             colData = coldata8,
                             design = ~sex)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("sex_M_vs_F")))
CMF_Colon <- subset(res, pvalue < 0.05 & log2FoldChange > 0 )
CMF_Colon <- as.data.frame(CMF_Colon)
CFM_Colon <- subset(res, pvalue < 0.05 & log2FoldChange < 0 )
CFM_Colon <- as.data.frame(CFM_Colon)

write.table(CMF_Colon,"CMF_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(CFM_Colon,"CFM_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)


dim(ApcM_Colon)
dim(ApcF_Colon)
BP1<- as.data.frame(simplify(enrichGO(gene=rownames(ApcM_Colon),keyType = "ENSEMBL",
                                      OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                      pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)))

BP1<- enrichGO(gene=rownames(ApcM_Colon),keyType = "ENSEMBL",
                                      OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                      pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)
View(as.data.frame(BP1))

BP2<- as.data.frame(simplify(enrichGO(gene=rownames(ApcF_Colon),keyType = "ENSEMBL",
                                      OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                      pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)))
View(as.data.frame(BP2))
BP2<-enrichGO(gene=rownames(ApcF_Colon),keyType = "ENSEMBL",
         OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
         pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(ApcM_Colon),
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")
KEGG1<- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",
                                      organism = 'mmu',
                                      pvalueCutoff = 1,qvalueCutoff = 1)%>%as.data.frame()
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(ApcF_Colon),
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")
KEGG2<- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",
                                    organism = 'mmu',
                                    pvalueCutoff = 1,qvalueCutoff = 1)%>%as.data.frame()
View(KEGG2)

############################### C vs B  Colon ########################
colnames(BR_newcount)
Colon <- BR_newcount[,grep("Colon",colnames(BR_newcount))]
colnames(Colon)
data <- Colon[,c(1:6,7:12)]
colnames(data)
coldata <- strsplit2(colnames(data),"_")%>%as.data.frame(.)
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition","sex","tissue","replicate")
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

res <- merge(dd,res,all.x=T,by="ENSEMBL")
write.table(res,"res_CvsB.txt", sep = '\t', col.names = NA, quote = FALSE)

View(res)
C_B_Colon <- subset(res, pvalue < 0.05 & log2FoldChange > 0 )						
C_B_Colon <- as.data.frame(C_B_Colon)
C_B_Colon$Class <- "CRC"
B_C_Colon <- subset(res, pvalue < 0.05 & log2FoldChange < 0 )						
B_C_Colon <- as.data.frame(B_C_Colon)
B_C_Colon$Class <- "BR"
CRC_BR_Diff <- rbind(C_B_Colon,B_C_Colon)
View(CRC_BR_Diff)
write.table(CRC_BR_Diff,"CRC_BR_Diff.txt", sep = '\t', col.names = NA, quote = FALSE)
B_anti <- data.frame(gene=intersect(Apc_Colon$SYMBOL,C_B_Colon$SYMBOL))
B_pro <- data.frame(gene=intersect(WT_Colon$SYMBOL,B_C_Colon$SYMBOL))
write.table(B_anti,"B_anti.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(B_pro,"B_pro.txt", sep = '\t', col.names = NA, quote = FALSE)

write.table(C_B_Colon,"C_B_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(B_C_Colon,"B_C_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)

BP1<- as.data.frame(clusterProfiler::simplify(enrichGO(gene=C_B_Colon$ENSEMBL,keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                     pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)))
View(B1)
BP2<- as.data.frame(clusterProfiler::simplify(enrichGO(gene=B_C_Colon$ENSEMBL,keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                     pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)))

View(BP2)

library(clusterProfiler)
CC1<-clusterProfiler::simplify(enrichGO(gene=C_B_Colon$ENSEMBL,keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="CC",pAdjustMethod="BH",
                                     pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T))%>%as.data.frame()
View(CC1)
CC2<- clusterProfiler::simplify(enrichGO(gene=B_C_Colon$ENSEMBL,keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="CC",pAdjustMethod="BH",
                                     pvalueCutoff=1,qvalueCutoff=1,readable=T))%>%as.data.frame()

View(CC2)



MF1<-clusterProfiler::simplify(enrichGO(gene=C_B_Colon$ENSEMBL,keyType = "ENSEMBL",
                                        OrgDb= org.Mm.eg.db, ont="MF",pAdjustMethod="BH",
                                        pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T))%>%as.data.frame()
View(MF1)
MF2<- clusterProfiler::simplify(enrichGO(gene=B_C_Colon$ENSEMBL,keyType = "ENSEMBL",
                                         OrgDb= org.Mm.eg.db, ont="MF",pAdjustMethod="BH",
                                         pvalueCutoff=1,qvalueCutoff=1,readable=T))%>%as.data.frame()

View(MF2)

meta <- res[res$SYMBOL%in%c("Slc11a2","Lox","Alas2","Icos","Ctla4",
                            "S100a9","Tigit","Lrrc32","Il17a","Lcn2",
                            "Cxcl13","Cxcl12","Tnfsf14","Efnb2","Gabrg2","Cbln2",
                            "Nptx1","Add2","Adgrb3","Nrxn2"),]
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
                     label = c('up'='up (num=1236)','noSig'='noSig (num=33896)','down'='down (num=1412)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('Apc vs BR')+theme_classic()
p + geom_text_repel(data = meta,aes(x = log2FoldChange,y = -log10(pvalue),label = meta$SYMBOL),
                    force=20,color="black",size=6,point.padding = 0.5,hjust = 0.5,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                    segment.color="black",segment.size=0.2,nudge_y=1)+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))

colon_tpm <- tpm[,grep("Colon",colnames(tpm))]%>%as.data.frame()
head(colon_tpm)
colon_tpm$ENSEMBL <- rownames(colon_tpm)

colon_tpm <- merge(colon_tpm,sym,all.x=T,by="ENSEMBL")
colon_tpm <- na.omit(colon_tpm[!duplicated(colon_tpm$SYMBOL),])
rownames(colon_tpm) <- colon_tpm$SYMBOL
#colon_tpm$SYMBOL <- NULL
colon_tpm$ENTREZID <- NULL
colon_tpm$ENSEMBL <- NULL
pheatmap::pheatmap(cor(colon_tpm))
ggdata <- melt(colon_tpm)
ggdata$group <- gsub("B","BR",strsplit2(ggdata$variable,"_")[,1])%>%
  gsub("C","Apc",.)
ggdata$group <- factor(ggdata$group,levels = c("WT","Apc","BR"))
data <- ggdata[ggdata$SYMBOL=="",]
ggplot(data,aes(x = group,y = log2(value),fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#F2B4A9","#88ADA6","#C2D6E4"))+
  geom_signif(comparisons = list(c("WT", "Apc"),c("Apc","BR"),c("WT","BR")),
              map_signif_level = F,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative expression log2(TPM)")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(legend.position = 'none')+ggtitle("Hoxc11")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))



View(ggdata)
data <- ggdata[ggdata$SYMBOL=="Il16",]
ggplot(data,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#F2B4A9","#88ADA6","#C2D6E4"))+
  geom_signif(comparisons = list(c("WT", "Apc"),c("Apc","BR"),c("WT","BR")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative expression log2(TPM)")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(legend.position = 'none')+ggtitle("Il16")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))


#dim(WT_Colon)						
View(Apc_Colon)						
library(clusterProfiler)						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(C_B_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_CvsB_Colon <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                               organism = 'mmu',						
                    
                                                                          pvalueCutoff = 0.05)%>%as.data.frame()						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(B_C_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_BvsC_Colon <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                               organism = 'mmu',						
                                               pvalueCutoff = 0.05)%>%as.data.frame()						
					
############################# C vs B Male Colon ##########################
colnames(Colon)
data <- Colon[,c(4:6,10:12)]
colnames(data)
coldata <- strsplit2(colnames(data),"_")%>%as.data.frame(.)
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data,
                             colData = coldata,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_C_vs_B")))						
C_B_M_Colon <- subset(res, pvalue < 0.05 & log2FoldChange > 0 )						
C_B_M_Colon <- as.data.frame(C_B_M_Colon)						
B_C_M_Colon <- subset(res, pvalue < 0.05 & log2FoldChange < 0 )						
B_C_M_Colon <- as.data.frame(B_C_M_Colon)						
write.table(C_B_M_Colon,"C_B_M_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(B_C_M_Colon,"B_C_M_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)


#dim(WT_Colon)						
#View(Apc_Colon)						
library(clusterProfiler)						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(C_B_M_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_CvsB_Colon_M <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                               organism = 'mmu',						
                                               pvalueCutoff = 0.05)%>%as.data.frame()						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(B_C_M_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_BvsC_Colon_M <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                               organism = 'mmu',						
                                               pvalueCutoff = 0.05)%>%as.data.frame()						

############################# C vs B Female Colon ##########################
colnames(Colon)
data <- Colon[,c(1:3,7:9)]
colnames(data)
coldata <- strsplit2(colnames(data),"_")%>%as.data.frame(.)
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data,
                             colData = coldata,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_C_vs_B")))						
C_B_F_Colon <- subset(res, pvalue < 0.05 & log2FoldChange > 0 )						
C_B_F_Colon <- as.data.frame(C_B_F_Colon)						
B_C_F_Colon <- subset(res, pvalue < 0.05 & log2FoldChange < 0 )						
B_C_F_Colon <- as.data.frame(B_C_F_Colon)						
write.table(C_B_F_Colon,"C_B_F_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(B_C_F_Colon,"B_C_F_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)

library(clusterProfiler)						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(C_B_F_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_CvsB_Colon_F <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                                 organism = 'mmu',						
                                                 pvalueCutoff = 0.05)%>%as.data.frame()						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(B_C_F_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_BvsC_Colon_F <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                                 organism = 'mmu',						
                                                 pvalueCutoff = 0.05)%>%as.data.frame()						



############################# B M vs F Colon ##########################
colnames(Colon)
data <- Colon[,c(1:6)]
colnames(data)
coldata <- strsplit2(colnames(data),"_")%>%as.data.frame(.)
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data,
                             colData = coldata,
                             design = ~sex)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("sex_M_vs_F")))				
BMF_Colon <- subset(res, pvalue < 0.05 & log2FoldChange > 0 )						
BMF_Colon <- as.data.frame(BMF_Colon)						
BFM_Colon <- subset(res, pvalue < 0.05 & log2FoldChange < 0 )						
BFM_Colon <- as.data.frame(BFM_Colon)
write.table(BMF_Colon,"BMF_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(BFM_Colon,"BFM_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)


#dim(WT_Colon)						
#View(Apc_Colon)						
library(clusterProfiler)						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(BMF_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_BM_Colon <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                                 organism = 'mmu',						
                                                 pvalueCutoff = 0.05)%>%as.data.frame()						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(BFM_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_BF_Colon <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                                 organism = 'mmu',						
                                                 pvalueCutoff = 0.05)%>%as.data.frame()						









############################# WT vs B Male Colon ##########################
colnames(Colon)
data <- Colon[,c(4:6,16:18)]
colnames(data)
coldata <- strsplit2(colnames(data),"_")%>%as.data.frame(.)
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data,
                             colData = coldata,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_WT_vs_B")))						
WT_B_M_Colon <- subset(res, pvalue < 0.05 & log2FoldChange > 0 )						
wT_B_M_Colon <- as.data.frame(WT_B_M_Colon)						
B_WT_M_Colon <- subset(res, pvalue < 0.05 & log2FoldChange < 0 )						
B_WT_M_Colon <- as.data.frame(B_WT_M_Colon)						
#dim(WT_Colon)						
#View(ApWT_Colon)						
library(clusterProfiler)						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(WT_B_M_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_WTvsB_Colon_M <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                                 organism = 'mmu',						
                                                 pvalueCutoff = 0.05)%>%as.data.frame()						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(B_WT_M_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_BvsWT_Colon_M <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                                 organism = 'mmu',						
                                                 pvalueCutoff = 0.05)%>%as.data.frame()						

############################# WT vs B Female Colon ##########################
colnames(Colon)
data <- Colon[,c(1:3,13:15)]
colnames(data)
coldata <- strsplit2(colnames(data),"_")%>%as.data.frame(.)
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition","sex","tissue","replicate")
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data,
                             colData = coldata,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_WT_vs_B")))						
WT_B_F_Colon <- subset(res, pvalue < 0.05 & log2FoldChange > 0 )						
WT_B_F_Colon <- as.data.frame(WT_B_F_Colon)						
B_WT_F_Colon <- subset(res, pvalue < 0.05 & log2FoldChange < 0 )						
B_WT_F_Colon <- as.data.frame(B_WT_F_Colon)						

library(clusterProfiler)						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(WT_B_F_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_WTvsB_Colon_F <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                                 organism = 'mmu',						
                                                 pvalueCutoff = 0.05)%>%as.data.frame()						
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(B_WT_F_Colon),				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						
KEGG_BvsWT_Colon_F <- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",				
                                                 organism = 'mmu',						
                                                 pvalueCutoff = 0.05)%>%as.data.frame()						

select <- read.delim("F:/A_MetaboAnalysis/Black_rice/weeks_22/KEGG_select.txt")

select$RichFactor=as.numeric(select$Count)/ as.numeric(sub("/\\d+", "",select$BgRatio))						

select$group <- factor(select$group,levels = c("WT_C","C_WT","B_C","C_B"))
ggplot(select,aes(group,Description)) + geom_point(aes(color=pvalue,size=RichFactor))+						
  scale_colour_gradient(low="#FF603F",high= "#6997ED")+						
  labs(size="RichFactor",x="",y="",title="")+					
  theme(axis.ticks = element_blank())+
  scale_y_discrete(limits=rev(unique(select$Description)))+theme_classic()+	
  theme(axis.title =element_text(size = 20,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 20,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.title  = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))
################################# AOM RNA ########################
setwd("F:/A_MetaboAnalysis/AOM_RNA/")
############################### C vs B  Colon ########################
AOM_new_count <- read.delim("F:/A_MetaboAnalysis/AOM_RNA/AOM_new_count.txt", row.names=1, comment.char="#")
library(rtracklayer)
library(tidyverse)
gtf = rtracklayer::import("F:/A_MetaboAnalysis/Black_rice/weeks_22/Mus_musculus.GRCm39.104.gtf")
gtf = as.data.frame(gtf)
tra = gtf[gtf$type=="transcript",
          c("start","end","gene_id")]
glt = mutate(tra, length = end - start) %>%
  .[order(.$length),] %>% 
  filter(!duplicated(gene_id)) 
final_gene <- as.data.frame(glt[,3:4])
rownames(final_gene) <- final_gene$gene_id
final_gene <- final_gene[rownames(AOM_new_count),]
kb <- final_gene$length / 1000
rpk <- AOM_new_count/ kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)%>%as.data.frame()
head(tpm)
colnames(tpm) <- strsplit2(colnames(tpm),"[.]")[,7]
tpm$ENSEMBL <- rownames(tpm)
dd <- tpm
dd$ENSEMBL <- rownames(dd)


colnames(AOM_new_count) <- strsplit2(colnames(AOM_new_count),"[.]")[,7]
group <- strsplit2(colnames(AOM_new_count),"_")[,1] 
group
colnames(Colon)
data <- AOM_new_count
colnames(data)
coldata <- strsplit2(colnames(data),"_")%>%as.data.frame(.)
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition","sex","replicate")
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
res <- merge(dd,res,all.x=T,by="ENSEMBL")
write.table(res,"res_AOM_CvsB.txt", sep = '\t', col.names = NA, quote = FALSE)
View(res)

C_B_Colon <- subset(res, pvalue < 0.05 & log2FoldChange > 0 )						
C_B_Colon <- as.data.frame(C_B_Colon)
C_B_Colon$Class <- "AOM/DSS"
B_C_Colon <- subset(res, pvalue < 0.05 & log2FoldChange < 0 )						
B_C_Colon <- as.data.frame(B_C_Colon)
B_C_Colon$Class <- "BR"
AOM_BR_Diff <- rbind(C_B_Colon,B_C_Colon)
KEGG_C_Colon <- clusterProfiler::enrichKEGG(gene =C_B_Colon$ENTREZID,keyType = "ncbi-geneid",				
                                             organism = 'mmu',						
                                             pvalueCutoff = 0.05)%>%as.data.frame()						
KEGG_B_Colon <- clusterProfiler::enrichKEGG(gene =B_C_Colon$ENTREZID,keyType = "ncbi-geneid",				
                                             organism = 'mmu',						
                                             pvalueCutoff = 0.05)%>%as.data.frame()						

View(KEGG_C_Colon)
View(KEGG_B_Colon)
write.table(KEGG_C_Colon,"KEGG_C_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(KEGG_B_Colon,"KEGG_B_Colon.txt", sep = '\t', col.names = NA, quote = FALSE)

BP_C<- as.data.frame(simplify(enrichGO(gene=C_B_Colon$ENSEMBL,keyType = "ENSEMBL",
                                     OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                     pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)))
BP_B<- as.data.frame(simplify(enrichGO(gene=B_C_Colon$ENSEMBL,keyType = "ENSEMBL",
                                       OrgDb= org.Mm.eg.db, ont="BP",pAdjustMethod="BH",
                                       pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)))
View(BP_B)
View(BP_C)

meta <- res[res$SYMBOL%in%c("Comt","Prkn","Pnkd","Sult1a1","Maob","Agtr1a","Cyp2d22",
                            "Folr1","Pnpo","Abcd4","Btd","Gsto2","Vnn1",
                            "Cd79a","Il27ra","Fcer2a","Sh2d1a","Tnfrsf13b",
                            "Slamf6","Slamf1","Irf4","Il18"),]
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
                     label = c('up'='up (num=452)','noSig'='noSig (num=36394)','down'='down (num=869)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('AOM/DSS vs BR')+theme_classic()
p + geom_text_repel(data = meta,aes(x = log2FoldChange,y = -log10(pvalue),label = meta$SYMBOL),
                    force=20,color="black",size=6,point.padding = 0.5,hjust = 0.5,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                    segment.color="black",segment.size=0.2,nudge_y=1)+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))



tpm$ENSEMBL <- rownames(tpm)
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=tpm$ENSEMBL,				
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")						

tpm <- merge(tpm,sym,all.x=T,by="ENSEMBL")
tpm <- na.omit(tpm[!duplicated(tpm$SYMBOL),])
rownames(tpm) <- tpm$SYMBOL

#colon_tpm$SYMBOL <- NULL
tpm$ENTREZID <- NULL
tpm$ENSEMBL <- NULL
#pheatmap::pheatmap(cor(colon_tpm))
ggdata <- melt(tpm)
ggdata$group <- gsub("B","BR",strsplit2(ggdata$variable,"_")[,1])%>%
  gsub("C","AOM/DSS",.)
ggdata$group <- factor(ggdata$group,levels = c("AOM/DSS","BR"))

data <- ggdata[ggdata$SYMBOL%in%c("Slc39a14"),]
ggplot(data,aes(x = group,y = log2(value+1),fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#F2B4A9","#88ADA6","#C2D6E4"))+
  geom_signif(comparisons = list(c("AOM/DSS","BR")),
              map_signif_level = T,na.rm =T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative expression log2(TPM)")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(legend.position = 'none')+ggtitle("Tnfrsf8")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))

select <- read.delim("KEGG_select.txt")
select$RichFactor=as.numeric(select$Count)/ as.numeric(sub("/\\d+", "",select$BgRatio))						
select$group <- factor(select$Class,levels = c("AOM/DSS","BR"))
ggplot(select,aes(group,Description)) + geom_point(aes(color=pvalue,size=RichFactor))+						
  scale_colour_gradient(low="#FF603F",high= "#6997ED")+						
  labs(size="RichFactor",x="",y="",title="")+					
  theme(axis.ticks = element_blank())+
  scale_y_discrete(limits=rev(unique(select$Description)))+theme_classic()+	
  theme(axis.title =element_text(size = 20,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 20,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.title  = element_text(size=12,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))

 a=WT_Colon[,c(26,21)]
 b=B_C_Colon[,c(26,21)]
 c=merge(a,b,all.x=T,by="SYMBOL")
 colnames(c) <- c("SYMBOL","log2FC_WTvsApc","log2FC_BvsC")
 B_pro_table <- c[c$SYMBOL%in%B_pro$gene,]%>%na.omit(.)
 
 a=Apc_Colon[,c(26,21)]
 b=C_B_Colon[,c(26,21)]
 c=merge(a,b,all.x=T,by="SYMBOL")
 colnames(c) <- c("SYMBOL","log2FC_ApcvsWT","log2FC_CvsB")
 B_anti_table <- c[c$SYMBOL%in%B_anti$gene,]%>%na.omit(.)
 write.table(B_anti_table,"B_anti_table.txt", sep = '\t', col.names = NA, quote = FALSE)
 
 write.table(B_pro_table,"B_pro_table.txt", sep = '\t', col.names = NA, quote = FALSE)
