
setwd("F:/A_MetaboAnalysis/Black_rice/tumor_purity/")
#install.packages("estimate", repos="http://R-Forge.R-project.org")
library(estimate)
#BiocManager::install("biomaRt")
library("biomaRt")
dat <-read.delim("F:/A_MetaboAnalysis/Black_rice/weeks_22/colon_tpm.txt", row.names=1)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
MtoH <- getLDS(attributes = "mgi_symbol", # 要转换符号的属性，这里基因名（第3步是基因名）
               filters = "mgi_symbol", #参数过滤
               mart = mouse, #需要转换的基因名的种属来源，也就是第2步的mouse
               values = rownames(dat), #要转换的基因集
               attributesL = "hgnc_symbol", #要同源转换的目标属性，这里还是转为基因名，也可加其他
               martL = human, #要同源转换的目标种属，也就是第2步的human
               uniqueRows = TRUE)
head(MtoH)
head(dat)
dat$MGI.symbol <- rownames(dat)
dd <- merge(dat,MtoH,all.x=T,by="MGI.symbol")
dd <- na.omit(dd)
dd <- dd[!duplicated(dd$HGNC.symbol),]
rownames(dd) <- dd$HGNC.symbol
dd$MGI.symbol <- NULL
dd$HGNC.symbol <- NULL
dd
dat <- dd
estimate <- function(dat, pro){
  input.f <- paste0(pro, "_estimate_input.txt")
  output.f <- paste0(pro, "_estimate_gene.gct")
  output.ds <- paste0(pro, "_estimate_score.gct")
  write.table(dat, file=input.f, sep="\t", quote=F)
  library(estimate)
  filterCommonGenes(input.f = input.f,
                    output.f = output.f,
                    id = "GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds = output.ds,
                platform = "illumina")
  scores <- read.table(output.ds, skip = 2, header = T)
  rownames(scores) = scores[,1]
  scores <- t(scores[,3:ncol(scores)])
  return(scores)
}
pro <- "tumor"
tumor_scores <- as.data.frame(estimate(dat, pro))
tumor_scores$tumor_purity <- cos(0.6049872018 + 0.0001467884 * tumor_scores$ESTIMATEScore)

tumor_scores$group <- strsplit2(rownames(tumor_scores),"_")[,1]
tumor_scores$group <- factor(tumor_scores$group,levels = c("WT","C","B"))

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}
data1 <- tumor_scores %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(ImmuneScore)) %>% 
  ungroup()

ggplot(tumor_scores, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_npg() +
  labs(x = "", 
  y = 'ESTIMATE Immune Score') +
  geom_signif(comparisons = list(c("WT", "B"),c("WT","C"),c("B","C")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2) +
  theme_classic(base_size = 16)+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "right",legend.title = NULL)

ggplot(data1, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_npg() +
  labs(x = "", 
       y = 'ESTIMATE Immune Score') +
  geom_signif(comparisons = list(c("WT", "B"),c("WT","C"),c("B","C")),
              map_signif_level = F,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2) +
  theme_classic(base_size = 16)+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "right",legend.title = NULL)

data1 <- tumor_scores %>% 
  dplyr::group_by(group) %>% 
  mutate(value = remove_outliers(StromalScore)) %>% 
  ungroup()
ggplot(data1, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_npg() +
  labs(x = "", 
       y = 'ESTIMATE Stromal Score') +
  geom_signif(comparisons = list(c("WT", "B"),c("WT","C"),c("B","C")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2) +
  theme_classic(base_size = 16)+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "right",legend.title = NULL)

ImmuneScore <- tumor_scores$ImmuneScore 
names(ImmuneScore) <- colnames(dat)
ImmuneScore <- as.data.frame(ImmuneScore)%>%t()%>%as.data.frame()

library(psych)
library(reshape2)
dim(ImmuneScore)
dim(as.data.frame(dat))
cor <-corr.test(t(ImmuneScore),t(dat),method = "spearman",adjust= "none")
dim(cor)
cmt <-t(cor$r)%>%as.data.frame()
head(cmt)
colnames(cmt) <- "R"
cmt$pvalue <- t(cor$p)
View(cmt)
cmt_select <- subset(cmt,pvalue<0.05)
cmt_select <- subset(cmt_select,abs(R)>0.6)
View(cmt_select)
dd <- data.frame(Correlation=cmt_select$R,row.names = rownames(cmt_select))
dd
pheatmap::pheatmap(dd,cluster_cols = F,show_rownames = F,fontsize = 18,angle_col = 0)

BP_up<- clusterProfiler::simplify(enrichGO(gene=rownames(cmt_select[cmt_select$R>0,]),keyType = "SYMBOL",
                                      OrgDb= org.Hs.eg.db, ont="BP",pAdjustMethod="BH",
                                      pvalueCutoff=0.05,qvalueCutoff=0.05))
BP_down<- clusterProfiler::simplify(enrichGO(gene=rownames(cmt_select[cmt_select$R<0,]),keyType = "SYMBOL",
                                        OrgDb= org.Hs.eg.db, ont="BP",pAdjustMethod="BH",
                                        pvalueCutoff=0.05,qvalueCutoff=0.05))

View(BP_up@result)
View(BP_down@result)
?dotplot
dotplot(BP_up,showCategory = 10,label_format = 100)
dotplot(BP_down,showCategory = 10,label_format = 100)

# Update the clusterprofilier to the latest Github version ( the lastest version is 4.7.1.3)
devtools::install_local("D:/DOSE_3.26.1.tar.gz")
devtools::install_local("D:/GOSemSim-master.zip")
remotes::install_github("YuLab-SMU/clusterProfiler") 

# Establish a local KEGG database

# install the packages

devtools::install_local("D:/createKEGGdb-master.zip")
# import the library and create a KEGG database locally 
library(createKEGGdb)

species <-c("ath","hsa","mmu", "rno","dre","dme","cel")
createKEGGdb::create_kegg_db(species)
# You will get KEGG.db_1.0.tar.gz file in your working directory
# install the KEGG.db and import it
install.packages("KEGG.db_1.0.tar.gz", repos=NULL,type="source")
library(KEGG.db)
# add use_internal_data=T in your enrichKEGG function


sym<- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(cmt_select[cmt_select$R>0,]),
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_up<- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",
                                    organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                    use_internal_data=T)
View(KEGG_up@result)
sym<- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(cmt_select[cmt_select$R<0,]),
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_down<- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",
                                      organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                      use_internal_data=T)
View(KEGG_down@result)
write.table(KEGG_up@result, file="KEGG_up.txt", sep="\t", quote=F)

write.table(KEGG_down@result, file="KEGG_down.txt", sep="\t", quote=F)
KEGG_up_select <- read.delim("F:/A_MetaboAnalysis/Black_rice/tumor_purity/KEGG_up_select.txt", 
                             row.names=1)
#View(KEGG_up_select[KEGG_up_select$select=="1",])
KEGG_up_select <- KEGG_up_select[KEGG_up_select$select==1,]%>%na.omit(.)

KEGG_up_select$BgRatio
install.packages("aPEAR")
library(aPEAR)
set.seed(348934)
?enrichmentNetwork
enrichmentNetwork(KEGG_up_select,colorBy ='pvalue',colorType='pval',nodeSize = "Count")




y=new("enrichResult",
      result=KEGG_up_select)
y@result
?dotplot

dotplot(y,showCategory = 20,label_format = 100,color = "pvalue")+theme_classic()+
  theme(axis.title =element_text(size = 20, color = 'black'),
        axis.text.x =element_text(size = 20, color = 'black'))+	
  theme(axis.text.y = element_text(size = 20, color = 'black'),
        legend.text = element_text(size=14, color = 'black'),
        legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8, color = 'black'),
        axis.line.y=element_line(size=0.8, color = 'black'))


StromalScore <- tumor_scores$StromalScore
names(StromalScore) <- colnames(dat)
StromalScore <- as.data.frame(StromalScore)%>%t()%>%as.data.frame()

library(psych)
library(reshape2)
dim(ImmuneScore)
dim(as.data.frame(dat))
cor <-corr.test(t(StromalScore),t(dat),method = "spearman",adjust= "none")
dim(cor)
cmt <-t(cor$r)%>%as.data.frame()
head(cmt)
colnames(cmt) <- "R"
cmt$pvalue <- t(cor$p)
View(cmt)
cmt_select <- subset(cmt,pvalue<0.05)
cmt_select <- subset(cmt_select,abs(R)>0.6)
View(cmt_select)
dd <- data.frame(Correlation=cmt_select$R,row.names = rownames(cmt_select))
dd
pheatmap::pheatmap(dd,cluster_cols = F,show_rownames = F,fontsize = 18,angle_col = 0)

sym<- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(cmt_select[cmt_select$R>0,]),
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_up<- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",
                                      organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                      use_internal_data=T)
View(KEGG_up@result)
sym<- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(cmt_select[cmt_select$R<0,]),
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_down<- clusterProfiler::enrichKEGG(gene =sym$ENTREZID,keyType = "ncbi-geneid",
                                        organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                        use_internal_data=T)
View(KEGG_down@result)
write.table(KEGG_up@result, file="KEGG_up.txt", sep="\t", quote=F)

write.table(KEGG_down@result, file="KEGG_down.txt", sep="\t", quote=F)





####################################################################
colon_tpm$genes <- rownames(colon_tpm)
ggdata <- melt(colon_tpm)
library(limma)
ggdata$group <- strsplit2(ggdata$variable,"_")[,1]
ggdata$group <- factor(ggdata$group,levels = c("WT","C","B"))

#View(ggdata)
data <- ggdata[ggdata$genes=="Ctla4",]
ggplot(data,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.1,lwd=1)+
  geom_boxplot(lwd=1,fatten = 1,width=0.5)+theme_classic()+
  scale_fill_manual(values=c("#83BFE3","#F2B4A9","#88ADA6"))+
  geom_signif(comparisons = list(c("WT", "C"),c("C","B")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative expression (TPM)") +xlab("")+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(size = 20),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Ctla4")+ 
  theme(plot.title = element_text(size = 20, face = "italic"))

data <- ggdata[ggdata$genes=="Lag3",]
ggplot(data[-c(6,1),],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.1,lwd=1)+
  geom_boxplot(lwd=1,fatten = 1,width=0.5)+theme_classic()+
  scale_fill_manual(values=c("#83BFE3","#F2B4A9","#88ADA6"))+ylab("Relative expression (TPM)") +xlab("")+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(size = 20),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Lag3")+ 
  theme(plot.title = element_text(size = 20, face = "italic"))

data <- ggdata[ggdata$genes=="Pdcd1",]
ggplot(data[-c(6,1),],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.1,lwd=1)+
  geom_boxplot(lwd=1,fatten = 1,width=0.5)+theme_classic()+
  scale_fill_manual(values=c("#83BFE3","#F2B4A9","#88ADA6"))+
  geom_signif(comparisons = list(c("WT", "C"),c("C","B")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative expression (TPM)") +xlab("")+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(size = 20),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Pdcd1")+ 
  theme(plot.title = element_text(size = 20, face = "italic"))


data <- ggdata[ggdata$genes=="Tigit",]
ggplot(data,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.1,lwd=1)+
  geom_boxplot(lwd=1,fatten = 1,width=0.5)+theme_classic()+
  scale_fill_manual(values=c("#83BFE3","#F2B4A9","#88ADA6"))+
  geom_signif(comparisons = list(c("WT", "C"),c("C","B")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative expression (TPM)") +xlab("")+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(size = 20),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Tigit")+ 
  theme(plot.title = element_text(size = 20, face = "italic"))

data <- ggdata[ggdata$genes=="Vista",]
ggplot(data,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.1,lwd=1)+
  geom_boxplot(lwd=1,fatten = 1,width=0.5)+theme_classic()+
  scale_fill_manual(values=c("#83BFE3","#F2B4A9","#88ADA6"))+
  geom_signif(comparisons = list(c("WT", "C"),c("C","B")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative expression (TPM)") +xlab("")+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(size = 20),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Vista")+ 
  theme(plot.title = element_text(size = 20, face = "italic"))


n=t(scale(t(dat[rownames(dat)%in%rownames(cmt_select[cmt_select$R>0,]),])))
n[n>2]=2
n[n<-2]=-2
pheatmap::pheatmap(na.omit(n),cluster_cols = F)

dat <- read.delim("F:/A_MetaboAnalysis/AOM_RNA/tpm_AOM.txt", row.names=1)
MtoH <- getLDS(attributes = "mgi_symbol", # 要转换符号的属性，这里基因名（第3步是基因名）
               filters = "mgi_symbol", #参数过滤
               mart = mouse, #需要转换的基因名的种属来源，也就是第2步的mouse
               values = rownames(dat), #要转换的基因集
               attributesL = "hgnc_symbol", #要同源转换的目标属性，这里还是转为基因名，也可加其他
               martL = human, #要同源转换的目标种属，也就是第2步的human
               uniqueRows = TRUE)
head(MtoH)
head(dat)
dat$MGI.symbol <- rownames(dat)
dd <- merge(dat,MtoH,all.x=T,by="MGI.symbol")
dd <- na.omit(dd)
dd <- dd[!duplicated(dd$HGNC.symbol),]
rownames(dd) <- dd$HGNC.symbol
dd$MGI.symbol <- NULL
dd$HGNC.symbol <- NULL
dd
dat <- dd
estimate <- function(dat, pro){
  input.f <- paste0(pro, "_estimate_input.txt")
  output.f <- paste0(pro, "_estimate_gene.gct")
  output.ds <- paste0(pro, "_estimate_score.gct")
  write.table(dat, file=input.f, sep="\t", quote=F)
  library(estimate)
  filterCommonGenes(input.f = input.f,
                    output.f = output.f,
                    id = "GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds = output.ds,
                platform = "illumina")
  scores <- read.table(output.ds, skip = 2, header = T)
  rownames(scores) = scores[,1]
  scores <- t(scores[,3:ncol(scores)])
  return(scores)
}
pro <- "tumor"
tumor2_scores <- as.data.frame(estimate(dat, pro))
tumor2_scores$tumor_purity <- cos(0.6049872018 + 0.0001467884 * tumor2_scores$ESTIMATEScore)

############### single cell ####################################
cor_result = cor_result[cor_result$gene %in% rownames(seu@assays[["RNA"]]@counts),]
tumor_gene = cor_result$gene[cor_result$cor > 0.3]#0.25 857  0.3 356
calc_cbscore <- function(){
  c_score <- apply(seu@assays[["RNA"]]@counts[tumor_gene,], 2, function(x){
    table(x>0)["TRUE"]/length(tumor_gene)
  })
  e_score = apply(seu@assays[["RNA"]]@counts,2,function(x){
    sum(x[tumor_gene])/sum(x)
  })
  combind_score = c_score * e_score
  kk = 1/(-log(combind_score))
  return(kk)
}
seu$tumor_score = calc_cbscore()
p2<-FeaturePlot(object = seu,
                reduction = "tsne",
                features = "tumor_score",
                pt.size=0.8)+
  scale_colour_gradientn(colours=c(c('#302e82','#110AEA','#110AEA','white',"red","red","red","#b20000")),
                         limits = c(0, 0.4),
                         breaks = c(0, 0.1, 0.2, 0.3, 0.4),
                         labels = c(0, 0.1, 0.2, 0.3, 0.4)) +  ## 停在0.175处
  #ggtitle(label = "cluster1 blueprint score")+
  theme(axis.title=element_text(size=40,face="bold"),    ## x,y坐标轴标题字体
        axis.text=element_text(vjust=1,size=40,face = "bold"),    ## x,y坐标轴刻度线字体
        # plot.title = element_text(size=40),
        plot.title=element_blank(),
        legend.text=element_text(size=35),    ## 图例字体
        legend.title = element_blank(),
        legend.key.size = unit(1.4,"line"))



#################################### Cecum ###########################################
dat <-cecum_tpm

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
MtoH <- getLDS(attributes = "mgi_symbol", # 要转换符号的属性，这里基因名（第3步是基因名）
               filters = "mgi_symbol", #参数过滤
               mart = mouse, #需要转换的基因名的种属来源，也就是第2步的mouse
               values = rownames(dat), #要转换的基因集
               attributesL = "hgnc_symbol", #要同源转换的目标属性，这里还是转为基因名，也可加其他
               martL = human, #要同源转换的目标种属，也就是第2步的human
               uniqueRows = TRUE)
head(MtoH)
head(dat)
dat$MGI.symbol <- rownames(dat)
dd <- merge(dat,MtoH,all.x=T,by="MGI.symbol")
dd <- na.omit(dd)
dd <- dd[!duplicated(dd$HGNC.symbol),]
rownames(dd) <- dd$HGNC.symbol
dd$MGI.symbol <- NULL
dd$HGNC.symbol <- NULL
dd
dat <- dd
estimate <- function(dat, pro){
  input.f <- paste0(pro, "_estimate_input.txt")
  output.f <- paste0(pro, "_estimate_gene.gct")
  output.ds <- paste0(pro, "_estimate_score.gct")
  write.table(dat, file=input.f, sep="\t", quote=F)
  library(estimate)
  filterCommonGenes(input.f = input.f,
                    output.f = output.f,
                    id = "GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds = output.ds,
                platform = "illumina")
  scores <- read.table(output.ds, skip = 2, header = T)
  rownames(scores) = scores[,1]
  scores <- t(scores[,3:ncol(scores)])
  return(scores)
}
pro <- "tumor"
tumor_scores <- as.data.frame(estimate(dat, pro))
tumor_scores$tumor_purity <- cos(0.6049872018 + 0.0001467884 * tumor_scores$ESTIMATEScore)

tumor_scores$group <- strsplit2(rownames(tumor_scores),"_")[,1]
tumor_scores$group <- factor(tumor_scores$group,levels = c("WT","C","B"))
data1 <- tumor_scores %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(ImmuneScore)) %>% 
  ungroup()

ggplot(tumor_scores, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_npg() +
  labs(x = "", 
       y = 'ESTIMATE Immune Score') +
  geom_signif(comparisons = list(c("WT", "B"),c("WT","C"),c("B","C")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2) +
  theme_classic(base_size = 16)+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "right",legend.title = NULL)

##################################TCGA#####################################
dat <- read.delim("F:/A_MetaboAnalysis/public_CRC/tcga_tpm.txt")
estimate <- function(dat, pro){
  input.f <- paste0(pro, "_estimate_input.txt")
  output.f <- paste0(pro, "_estimate_gene.gct")
  output.ds <- paste0(pro, "_estimate_score.gct")
  write.table(dat, file=input.f, sep="\t", quote=F)
  library(estimate)
  filterCommonGenes(input.f = input.f,
                    output.f = output.f,
                    id = "GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds = output.ds,
                platform = "illumina")
  scores <- read.table(output.ds, skip = 2, header = T)
  rownames(scores) = scores[,1]
  scores <- t(scores[,3:ncol(scores)])
  return(scores)
}
pro <- "tumor"
tumor_scores <- as.data.frame(estimate(dat, pro))
tumor_scores$tumor_purity <- cos(0.6049872018 + 0.0001467884 * tumor_scores$ESTIMATEScore)
View(tumor_scores)
meta<- read.delim("F:/A_MetaboAnalysis/public_CRC/COAD_meta_data.txt")
meta <- COAD_meta_data[colnames(dat),]
tumor_scores$group <- meta$tumortype

table(meta$tumortype)





ggplot(tumor_scores, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_npg() +
  labs(x = "", 
       y = 'ESTIMATE Immune Score') +
  geom_signif(comparisons = list(c("Tumor","Normal")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2) +
  theme_classic(base_size = 16)+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "right",legend.title = NULL)
table(meta$CMS_type)

tumor_scores$CMS_type <- paste0("CMS",meta$CMS_type)
which(meta$tumortype=="Normal")
tumor_scores$CMS_type[1:41] <- "Normal"
tumor_scores <- tumor_scores[grep("CMSNA",tumor_scores$CMS_type,invert = T),]
tumor_scores$CMS_type <- factor(tumor_scores$CMS_type,levels = c("Normal","CMS1","CMS2","CMS3","CMS4"))
ggplot(tumor_scores, aes(x = CMS_type, y = ImmuneScore, fill = CMS_type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_npg() +
  labs(x = "", 
       y = 'ESTIMATE Immune Score') +
  geom_signif(comparisons = list(c("CMS2","CMS4")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2) +
  theme_classic(base_size = 16)+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "right",legend.title = NULL)

tumor_scores$tumor_purity
ggplot(tumor_scores, aes(x = CMS_type, y = tumor_purity, fill = CMS_type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_npg() +
  labs(x = "", 
       y = 'Estimate Tumor Purity Score') +
  geom_signif(comparisons = list(c("CMS2","CMS4")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2) +
  theme_classic(base_size = 16)+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "right",legend.title = NULL)

tumor_scores$StromalScore
ggplot(tumor_scores, aes(x = CMS_type, y = StromalScore, fill = CMS_type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_npg() +
  labs(x = "", 
       y = 'Estimate Stromal Score') +
  geom_signif(comparisons = list(c("CMS2","CMS4")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2) +
  theme_classic(base_size = 16)+	
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "right",legend.title = NULL)
