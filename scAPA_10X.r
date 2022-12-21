# NCBI下载10X单细胞测序数据
prefetch SRR84275{19..26}

#数据预处理python脚本
import os
data_dir="SRA"
SAMPLES=os.listdir(data_dir)
i=0
while i<len(SAMPLES):
        SAMPLES[i]=SAMPLES[i].split('.')[0]
        i+=1
samples=list(set(SAMPLES))
rule all:
    input:    
        expand("fastq/{sample}_S1_L001_I1_001.fastq.gz", sample=samples),
        expand("fastq/{sample}_S1_L001_R1_001.fastq.gz", sample=samples),
        expand("fastq/{sample}_S1_L001_R2_001.fastq.gz", sample=samples)
rule fastqdump:
    input:
         "SRA/{sample}.sra"
    output:
         "fastq/{sample}_1.fastq",
         "fastq/{sample}_2.fastq",
         "fastq/{sample}_3.fastq"
    log:
        "logs/{sample}.SRA2fastq.log"
    shell:
        "/home/wangby/.conda/envs/wby/bin/fastq-dump --split-files -O fastq/ {input} >{log} 2>&1"
rule pigz:
    input:
        "fastq/{sample}_1.fastq",
        "fastq/{sample}_2.fastq",
        "fastq/{sample}_3.fastq"
    output:
        "fastq/{sample}_1.fastq.gz",
        "fastq/{sample}_2.fastq.gz",
        "fastq/{sample}_3.fastq.gz"
    log:
        "logs/{sample}.fastq2fqgz.log"
    shell:
        "pigz -p 20 {input} >{log} 2>&1"
rule changenames:
    input:
        "fastq/{sample}_1.fastq.gz",
        "fastq/{sample}_2.fastq.gz",
        "fastq/{sample}_3.fastq.gz"
    output:
        "fastq/{sample}_S1_L001_I1_001.fastq.gz",
        "fastq/{sample}_S1_L001_R1_001.fastq.gz",
        "fastq/{sample}_S1_L001_R2_001.fastq.gz"
    shell:
        """
        mv {input[0]} {output[0]}
        mv {input[1]} {output[1]}
        mv {input[2]} {output[2]}
        """

#使用cellranger进行比对
ls ../SRA |while read i;do i=${i:0:10};echo $i;
/opt/cellranger-6.1.2/cellranger count \
--transcriptome=/DATA/public/refdata-gex-GRCh38-2020-A \
--fastqs=../fastq \
--id=$i \
--sample=$i \
--localcores=30 >../logs/${i}.cellranger.log;
mkdir -p ../expma_bam_html/$i/
cp -r $i/outs/filtered_feature_bc_matrix/* $i/outs/possorted_genome_bam.bam \
$i/outs/web_summary.html ../expma_bam_html/$i/

######使用scAPAtrap进行APA位点的探查######

# 说在前面 --------------------------------------------------------------------
####
# refined by wangby 2022-04-16 | wangboyuan@mail.nankai.edu.cn
# 因为scAPAtrap最慢的步骤都是单线程，所以能够一起并行跑还是尽量去并行跑
# 目前的解决方案，个人觉得很有效的是在shell界面放在后台一起跑，切到样本目录上执行如下命令
# ls |while read i;do cd $i;nohup Rscript 000scAPAtrap_upstream.R &; cd ../ done
# 该脚本仅仅适用于10X测序文件
####

# Step 0:准备APA ---------------------------------------------------------------
# 输入软件绝对路径
samtools.path <- '/home/wangby/.conda/envs/wby/bin/samtools'
umitools.path <- '/home/wangby/.conda/envs/wby/bin/umi_tools'
featureCounts.path <- "/home/wangby/.conda/envs/wby/bin/featureCounts"
star.path <- '/home/wangby/.conda/envs/wby/bin/STAR'
# 设置线程数
threads <- 10
# 设置物种 MMU/HSA
species="MMU"
# 设置最大的peak宽度，默认为1000，实际上一般peak的中位数在200-300左右
maxwidth <- 1000
# R2读长，由输入文件的bam中每条reads的读长决定，不同试剂盒和测序公司结果不一
readlength <- 98
# bam文件的RNAME列是否有chr前缀，默认TRUE
IsChr = T
# 输出路径，一般默认当前路径就好
outputdir <- '.'
# 输入bam文件路径
demo.bam <- './possorted_genome_bam.bam'
# 输入barcode压缩文件，默认为 ./barcodes.tsv.gz
barcode <- "./barcodes.tsv.gz"
# 找到的peak和tails被认为可以视作一个PAs的碱基距离，默认是50，可以设置为0
d <- 50
# 加载scAPAtrap包
library(scAPAtrap)
# bam文件的RNAME一列是否有chr前缀的正则表达式转换
# system("samtools view -H ./bamfolder/pbmc1k.bam|sed -e 
#        's/chr\([0-9XY]\)/\1/' -e 's/chrM/MT/' |samtools reheader 
#        - ./bamfolder/pbmc1k.bam > ./bamfolder/pbmc1k.removeCHR.bam")
# Step 1: findUniqueMap ---------------------------------------------------
# 去除FlAG为256的reads：即比对到参考基因组多处的序列，否则会对后续分析造成影响
# sort&&index 下一步umi_tools需要
nextinput <- findUniqueMap(samtools.path, input = demo.bam, thread = threads, index = T, sort = T)
# Step 2: dedupByPos ------------------------------------------------------
# 去除CB，UB，read完全一样的reads————保留质量最好的reads，较慢
nextinput <- dedupByPos(umitools.path, nextinput, TenX = T)
# Step 3: separateBamBystrand ---------------------------------------------
# 分开±链，FLAG为16是反向比对，0是正向
nextinput <- separateBamBystrand(samtools.path, nextinput, 10)
# Step 4:findPeaks --------------------------------------------------------
# 设置染色体索引（应该直接访问bam文件，用正则表达式提出来）
if (species =="HSA"||species == "MMU"){
  if (species == "HSA"){
    if (IsChr) {
      chrs <- paste0("chr",c(as.character(1:22),'X','Y'))
      
    }
    chrs <- c(as.character(1:22),'X','Y')
  }
  if (species == "MMU"){
    if (IsChr) {
      chrs <- paste0("chr",c(as.character(1:19),'X','Y'))
    }
    chrs <- c(as.character(1:19),'X','Y')

  }
}else {
        print("please input the right species:HSA or MMU")
}
# 首先，使用loadBpCoverages函数导入整个BAM文件，并计算染色体每个碱基的覆盖率。
fullcovF <- loadBpCoverages(nextinput[1],chrs)
fullcovR <- loadBpCoverages(nextinput[2],chrs)
saveRDS(fullcovF,"fullcovF.rds")
saveRDS(fullcovR,"fullcovR.rds")
# 然后，使用FindPeak根据计算出的覆盖率识别峰值（较慢）
forwardPeaks <-findPeaks(fullcovF, '+', readlength, maxwidth)
reversePeaks <-findPeaks(fullcovR, '-', readlength, maxwidth)
saveRDS(forwardPeaks,"forwardPeaks.rds")
saveRDS(reversePeaks,"reversePeaks.rds")
# 最后，使用generateSAF生成一个SAF文件，记录已识别峰值的信息。
peaksfile <- generateSAF(forwardPeaks, reversePeaks, outputdir)
# Step 5:countPeaks -------------------------------------------------------
# 用上一步得到的peak文件注释原始bam文件
final.bam <- generateFinalBam(featureCounts.path,samtools.path,demo.bam,peaksfile,thread = threads)
# 使用注释过的final BAM文件计算每个峰值的读取数
counts.tsv <- countPeaks(umitools.path,final.bam,outputdir,TenX=T)
# Step 6:FindTails --------------------------------------------------------
#识别poly（A）站点所在的潜在区域。
#精确定位poly（A）站点
tails <- findTails(bamfile = demo.bam)
saveRDS(tails,"tails.rds")

# 使用Seurat对样本进行质控
library(Seurat)
library(dplyr)
files <- list.files("./data/000.rawX8/")
sscSeuratList <- list()
for (i in 1:length(files)) {
  expr <- Read10X(data.dir = paste0("./data/000.rawX8/",files[i]))
  seurat_obj <- CreateSeuratObject(
    counts = expr,
    project = files[i],
    min.cells = 3,
    min.features = 200 
  )  
  #质控和选择细胞
  seurat_obj[["percent.mt"]] <-
    PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj <- subset(seurat_obj, subset =
                         nFeature_RNA > 1000 &
                         nFeature_RNA < 6000 &
                         percent.mt < 25)
  seurat_obj <- NormalizeData(seurat_obj,verbose = FALSE) %>%
                FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  seurat_obj <- RunPCA(seurat_obj,features = VariableFeatures(object = seurat_obj),vars.to.regress = "percent.mt")
  seurat_obj@meta.data$group <- files[i]
  sscSeuratList[[i]] <- seurat_obj
  print(files[i])
}
saveRDS(sscSeuratList,"./data/004.intergrated_seurat/sscSeuratList.rds")

##==harmony整合多样本==##
library(harmony)
library(ggplot2)
scRNAlist <- readRDS("./data/004.intergrated_seurat/sscSeuratList.rds")
##PCA降维
scRNA_harmony <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], scRNAlist[[4]], scRNAlist[[5]], 
                                           scRNAlist[[6]], scRNAlist[[7]], scRNAlist[[8]]))
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
##整合
system.time({scRNA_harmony <- RunHarmony(seurat_object, group.by.vars = "orig.ident",lambda =2)})
DimPlot(seurat_object,reduction = "harmony", group.by = "orig.ident",label=T)
DimPlot(scRNA_harmony,reduction = "harmony", group.by = "orig.ident",label=T)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- FindClusters(scRNA_harmony,resolution = 0.6)
DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident",label = T,repel = T)+
DimPlot(scRNA_harmony,reduction = "harmony",label = T)
scRNA_harmony <- RunUMAP(scRNA_harmony,reduction = "harmony",dims = 1:30,resolution =0.6)
DimPlot(scRNA_harmony,reduction = "umap",label = T)+
  DimPlot(scRNA_harmony,reduction = "umap",group.by = "group",label = T)
scRNA_harmony <- RunTSNE(scRNA_harmony,reduction = "harmony",dims = 1:30)
DimPlot(scRNA_harmony,reduction = "tsne",label = T)+
  DimPlot(scRNA_harmony,reduction = "tsne",group.by = "group",label = T)
saveRDS(scRNA_harmony,"./data/007.finalPACds&seuratobject/Fialseurat_object.rds")
# 按需去除不需要的亚群
# 去除重复的barcodes并且写入到每个样本的文件中去
all_barcodes_rmdup <- data_frame()
seurat_object <- readRDS("./data/004.intergrated_seurat/sscSeurat_harmony_filteredbyDdx4.rds")
all_barcodes <- stringr::str_split(colnames(seurat_object),"_",simplify = TRUE)[(!duplicated(stringr::str_split(colnames(seurat_object),"_",simplify = TRUE)[,1])),1]
samples <- list.files("data/000.rawX8/")
for (i in 1:length(samples)) {
  print(i)
  barcode <- read.delim2(paste0("./data/000.rawX8/",samples[i],"/barcodes.tsv"),header = F)
  barcode <- barcode$V1[barcode$V1 %in% all_barcodes]
  write.table(barcode,paste0("./data/000.rawX8/",samples[i],"/barcode_rmdup.csv"))
}
# 对得到的APA矩阵进行注释

# 人
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene",force = TRUE)
# txdbhg38 <- parseGenomeAnnotation(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
# table(txdbhg38$anno.need$type)
# save(txdbhg38, file='txdbhg38.rda')

# 小鼠
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene",force = TRUE)
# txdbmm10 <- movAPA::parseGenomeAnnotation(TxDb.Mmusculus.UCSC.mm10.ensGene::TxDb.Mmusculus.UCSC.mm10.ensGene)
# table(txdbmm10$anno.need$type)
# saveRDS(txdbmm10, file='./data/gff.anno/txdbmm10.rds')

expmalist <- paste0("./data/001.expma_d=50/",list.files("./data/001.expma_d=50/"))
for (i in expmalist) {
  print(basename(i))
  expma <- readRDS(i)
  coldata <- data.frame(group = colnames(expma)[7:ncol(expma)], row.names = colnames(expma)[7:ncol(expma)])
  scPACds <- movAPA::readPACds(pacFile = expma, colDataFile = coldata)
  saveRDS(scPACds,paste0("./data/002.50_anno_scPACds/scPACds_",
                         stringr::str_split(list.files("./data/001.expma_d=50/"),"50_",simplify = TRUE)[,2][which(expmalist[]==i)]))
}
txdbmm10 <- readRDS("./data/000.gff.anno/txdbmm10.rds")

# 有几点需要注意一下，在执行该代码的时候，总是报错找不到bsgenome这个对象，原因未知，
# 最后直接start with local environment，并提前在环境里加载了bsgenome 解决
# 要注意的是，这一步非常吃内存，所以，尽量把IP去了，减少不必要的内存占用，否则700G的内存也无法满足需求

library(BSgenome)
library(movAPA)
library(BSgenome.Mmusculus.UCSC.mm10)
PACdsList <- list()
filelist <- list.files("./data/001.expma_d=50_BCuniq/")
bsgenome <- BSgenome.Mmusculus.UCSC.mm10
gff <- movAPA::parseGff("./data/000.gff.anno/gencode.vM28.annotation.gff3")
for (i in 1:length(filelist)) {
  print(filelist[i])
  expma <- readRDS(paste0("./data/001.expma_d=50_BCuniq/",filelist[i]))
  coldata <- data.frame(group = colnames(expma)[7:ncol(expma)], row.names = colnames(expma)[7:ncol(expma)])
  PACdsList[[i]] <- movAPA::readPACds(pacFile = expma, colDataFile = coldata)
  PACdsList[[i]] <- removePACdsIP(PACdsList[[i]], bsgenome, returnBoth = FALSE, 
                                  up=-10, dn=10, conA=6, sepA=7)
  PACdsList[[i]]@colData$sample <- stringr::str_split(filelist[i],"_",simplify = TRUE)[,2]
}
mergedPACds <- mergePACds(PACdsList, d=24)
mergedPACds@anno <- movAPA::annotatePAC(mergedPACds, txdbmm10)@anno[rownames(mergedPACds@counts),] 
#使mregedPACds中的细胞与Seurat对象中的保持一致
mergedPACds <- subsetPACds(mergePACds,pool = FALSE,group = "group",conds = allbarcodes)
saveRDS(mergedPACds,"./data/006.mergedPACds/mergedPACds.rds")
#处理得到的整合后的scPACds对象
library(movAPA)
# 先对得到的PACdataset进行质控
# 1.拓展3’UTR的PACs
scPACds=ext3UTRPACds(scPACds, ext3UTRlen=2000)
#让我们把目光集中在ext3‘UTR上吧。
# 1.去掉内含子和polyA Clusters<20的PAs
scPACdsFlt=subsetPACds(scPACds, totPACtag=20, choosePA=NULL,
                       noIntergenic=TRUE, verbose=TRUE) 
# 2'.只保留3'UTR的PACs，并且大于两个PACs
scPACdsFlt=get3UTRAPAds(scPACdsFlt, sortPA=TRUE, choose2PA=NULL)
# 2.数据标准化
scPACdsNorm2=movAPA::normalizePACds(scPACdsFlt, method='TPM') 
saveRDS(scPACdsNorm2,"./data/007.finalPACds&seuratobject/ext3nor_scPACds.rds")
# 为scPACdsNorm2添加补充信息
test <- scRNA_harmony@reductions$umap@cell.embeddings
test1 <- scRNA_harmony@reductions$tsne@cell.embeddings
test <- merge(test,test1,by = "row.names")
test1 <- scRNA_harmony@meta.data
test1$Row.names <-rownames(test1) 
test2 <- inner_join(test,test1,by = "Row.names")
test2 <- test2[,-c(7,8,9,10,12,13)]
test2$Row.names <- stringr::str_split(test2$Row.names,"-",simplify = T)[,1]
rownames(test2) <- test2$Row.names
test2$Row.names <- NULL
test2$group <- rownames(test2)
scPACdsNorm2@colData <- left_join(scPACdsNorm2@colData,test2,by = "group")

# 计算单细胞水平APA的统计信息
scPACdsStat=movStat(scPACds, minPAT=c(1,5), ofilePrefix='scPACds.stat')
scPACdsStat$pat5$group <- rownames(scPACdsStat$pat5)
inner_join(scPACdsStat$pat5,scPACdsNorm2@colData,by = "group")
ggplot(data = test3,aes(x=UMAP_1,y=UMAP_2,colour = log(nPAC)))+
geom_point()+
scale_color_gradient2(low="yellow",mid="red",high="black",midpoint = 7.5)
ggplot(data = test3,aes(x=UMAP_1,y=UMAP_2,colour = log(nGene)))+
geom_point()+
scale_color_gradient2(low="yellow",mid="red",high="black",midpoint = 7.5)

# 计算样本水平的APA统计信息
scPACdsNorm2SA <- subsetPACds(scPACdsNorm2,group = "sample",pool = TRUE)
pstatsSA=movStat(scPACdsNorm2SA, minPAT=c(5, 10), ofilePrefix=NULL)
plotPACdsStat(pstats, pdfFile='PACds_stat.pdf', minPAT=c(5,10))

#  计算总体PAS信号分布并可视化
scPACdsNorm2PAS = annotateByPAS(scPACdsNorm2, bsgenome, grams='V1',
from=-50, to=-1, label=NULL)
table(scPACdsNorm2PAS@anno$V1_gram)
pas=scPACdsNorm2PAS@anno$V1_gram[!is.na(scPACdsNorm2PAS@anno$V1_gram)]
plotSeqLogo(pas)
#  计算PAS上下游碱基分布情况

faFiles=faFromPACds(scPACdsNorm2, bsgenome, what='updn', fapre='updn',
up=-300, dn=100, byGrp='ftr')
faFiles=c("updn.3UTR.fa", "updn.CDS.fa", "updn.exon.fa", "updn.intron.fa")
plotATCGforFAfile (faFiles, ofreq=FALSE, opdf=FALSE,
refPos=301, mergePlots = TRUE)

#  计算单个样本的PAS信号分布并可视化
E16_5PACds <- subsetPACds(scPACdsNorm2,group = "sample",cond1 = "E16.5")
E16_5PACdsPAS=annotateByPAS(E16_5PACds, bsgenome, grams='V1',
from=-50, to=-1, label=NULL)
table(E16_5PACdsPAS@anno$V1_gram)
pas=E16_5PACdsPAS@anno$V1_gram[!is.na(E16_5PACdsPAS@anno$V1_gram)]
plotSeqLogo(pas)

#  绘制不同marker在不同细胞亚群的表达气泡图
DotPlot(seurat_object,features= c("Dppa3","Nanog","Sox2","Etv5","Id4","Lhx1","Ret","Dnmt3a","Sohlh1","Neurog3","Sox3","Kit","Stra8","Ddx4"))

# 绘制不同marker在不同细胞亚群中表达情况的UMAP降维图
FeaturePlot(seurat_object,features = c("Cpsf1","Cpsf2","Cpsf3","Cpsf4","Wdr33","Cstf1","Cpsf6","Cpsf7","Pcf11","Clp1","Rbbp6","Pabpn1"),ncol = 3)

#  寻找关键亚群之间的差异APA事件的基因并进行功能富集分析
PACds90 <- subsetPACds(scPACdsNorm2,group = "seurat_clusters",cond1 = "9",cond2 = "0",pool = FALSE)
PACds90 = mergedPACds(PACds90,d = 100)
sw90=movAPAswitch(PACds=PACds90 , group='seurat_clusters',
avgPACtag=0, avgGeneTag=0,
only3UTR=TRUE, mergeReps='pool',
aMovDEPACRes=NULL, DEPAC.padjThd=NULL, nDEPAC=0,
mindist=0, fisherThd=0.05, logFCThd=0, cross=FALSE,
selectOne='fisherPV')
heat90=movRes2heatmapResults(sw90)
heat90=subsetHeatmap(heat90, padjThd=0.001, valueThd=2)
heat90@value=heat90@value[rowSums(is.na(heat90@value))==0, ]
heat90@value=heat90@value[order(rowMeans(heat90@value), decreasing =T ), ]

#  加载go分析所用的R包
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(DOSE)

gene90 <- bitr(rownames(heat90), fromType = "ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb = "org.Mm.eg.db")
go90long <- enrichGO(gene90$ENTREZID[1:100],
                        keyType = 'ENTREZID',
                        OrgDb = org.Mm.eg.db, 
                        ont='ALL', 
                        pAdjustMethod = 'BH',
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.2,
                        readable = TRUE) 
                        showCategory = 15,
                        drop = T)

length90=length(gene90$ENTREZID)
go90long <- enrichGO(gene90$ENTREZID[(length90-99):length90],
                        keyType = 'ENTREZID',
                        OrgDb = org.Mm.eg.db, 
                        ont='ALL', 
                        pAdjustMethod = 'BH',
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.2,
                        readable = TRUE) 

barplot(go90short,
        showCategory = 15,
        drop = T)

barplot(go90long,
        showCategory = 15,
        drop = T)

