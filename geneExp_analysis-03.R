## Author: Kelly Quek
## Date: 2015-Dec-14
## Last update: 2016-Jan-06
## Objective: To analyse gene expression for ALL Teresa project (removed outliers)

# load libraries
library(Rsubread)
library(limma)
library(edgeR)

# read in target file
options(digits=2)
targets <- readTargets(file="sample_target-03.txt")

# count numbers of reads mapped to NCBI Refseq genes
fc <- featureCounts(files=targets$OutputFile_align,annot.inbuilt="hg19", isPairedEnd = TRUE, nthreads=24, ignoreDup=TRUE, reportReads=TRUE)

#count matrix
x <- DGEList(counts=fc$counts, genes=fc$annotation, group=targets$group)

# generate RPKM values if you need them
x_rpkm <- rpkm(x,x$genes$Length)

# filter out low-count genes. Only include genes where expression in cpm (count per million)
### is greater than 1 in at least 1 sample
isexpr <- rowSums(cpm(x) > 1) >= 1
x <- x[isexpr,]

y <- calcNormFactors(x)


## Gene annotation
ncbi.L1 <- readLines("/Users/kelly.quek/Downloads/references/Homo_sapiens.gene_info", n=1)
ncbi.colname <- unlist(strsplit(substring(ncbi.L1, 10, 234), ' '))
ncbi <- read.delim("/Users/kelly.quek/Downloads/references/Homo_sapiens.gene_info", skip=1, header=FALSE, stringsAsFactors=FALSE)
colnames(ncbi) <- ncbi.colname
m <- match(y$genes$GeneID, ncbi$GeneID)
y$genes$Chr <- ncbi$chromosome[m]
y$genes$Symbol <- ncbi$Symbol[m]
y$genes$Strand <- NULL
head(y$genes)


# create a design matrix
#type <- factor(targets$group)
type <- factor(targets$TLDA)
design <- model.matrix(~type)  # might need to change
design

# cluster libraries
plotMDS(y, col=c(type))

## voom with quality weight (I used this due to sample quality observed from MDS plot)
vwts <- voomWithQualityWeights(y, design=design, normalization="none", plot=TRUE)
vfit2 <- lmFit(vwts)
vfit2 <- eBayes(vfit2)

top2 <- topTable(vfit2,coef=2,sort.by="P", num=Inf)
sum(top2$adj.P.Val<0.05)

top3 <- topTable(vfit2,coef=3,sort.by="P", num=Inf)
sum(top3$adj.P.Val<0.05)

# Compute counts per million
# this step is to generate the counts with 0
z <- cpm(y)
z <- as.data.frame(z)
z$GeneID <- rownames(z)
cz <- cbind(z, y$counts)
label <- read.table("./label.txt", header = TRUE, sep = "\t")
colnames(cz) <- label$x 

#total.2 <- merge(top2, cz, by.x="GeneID", by.y="GeneID")
#write.table(total.2, file="Teresa03_Neg_vs_Pos_limma_voom_weight.csv", sep=",", row.names=FALSE)

total.2 <- merge(top2, cz, by.x="GeneID", by.y="GeneID")
write.table(total.2, file="Teresa03_grp1_vs_grp2_limma_voom_weight.csv", sep=",", row.names=FALSE)

total.3 <- merge(top3, cz, by.x="GeneID", by.y="GeneID")
write.table(total.3, file="Teresa03_gpr1_vs_grp3_limma_voom_weight.csv", sep=",", row.names=FALSE)

h1 <-  topTable(vfit2,coef=2,sort.by="P", num=Inf, p.value=0.05) 
h1 <-  topTable(vfit2,coef=3,sort.by="P", num=Inf, p.value=0.05) 
#h1 <-  topTable(vfit2,coef=2,sort.by="P", num=100, p.value=1)
#h2 <-  topTable(vfit2,coef=3,sort.by="P", num=100, p.value=1)


#### For the heatmap, the counts are log transformed.
a <- cpm(y, log=TRUE)
colnames(a) <- targets$patient_id
a1 <- a[which(rownames(a) %in% h1$GeneID),]

# set column names of a1 as group
colnames(a1) <- targets$group

library(gplots)
## heatmap
heatmap.2(a1, col=greenred(75), scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=0.8)

# save image
save.image(file = "teresa_GeneExp-03.RData")


######## Volcano Plot ##########
#
res <- read.csv("Teresa03_Neg_vs_Pos_limma_voom_weight.csv", sep = ",", header = TRUE)
head(res)
names(res)
min(res$log2FoldChange)
max(res$log2FoldChange)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(P.Value), pch=1, main="Volcano plot\nTeresa03_Neg_vs_Pos", xlim=c(-7.0,7.5)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
#with(subset(res, adj.P.Val<.05 ), points(log2FoldChange, -log10(P.Value), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(P.Value), pch=20, col="orange"))
#with(subset(res, adj.P.Val<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(P.Value), pch=20, col="green"))
with(subset(res, adj.P.Val<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(P.Value), pch=20, col="red"))
# Label points with the textxy function from the calibrate plot
library(calibrate)   # label a few genes of interest 
with(subset(res, adj.P.Val<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(P.Value), labs=Symbol, cex=.8))



