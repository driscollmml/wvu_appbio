setwd("C:/Users/arpol/Not_Bad/R_packages/")
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

#read in the folder names with the read data
pheno_data = read.csv("run_data.csv")

#here we specifically made the files into a ballgown variable
bg_chrX = ballgown(dataDir = "ballgown", samplePattern = "lane", pData=pheno_data)

#extract transcript with both coverage and fpkm
Whole_bg_chrx=texpr(bg_chrX,'all')

##Extract expression (expr) values (as FPKM) for  genes (g). 
gene_expression = gexpr(bg_chrX)

#extract expression for genes in fpkm and coverage
gene_bg_chrx=gexpr(bg_chrX)

##Look at the exon, intron and transcript  data
structure(bg_chrX)$exon

#store mapping between transcripts and genes
transcript_gene_table=indexes(bg_chrX)$t2g

#check the unique number of genes in the transcript/genes table
length(unique(transcript_gene_table[,"g_id"]))

#plot average transcript length
hist(Whole_bg_chrx$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col = "darkslategray1")

#how many transcripts are there per gene? count the number of genes and count the number of transcripts pere gene and plot it.
counts=table(transcript_gene_table[,"g_id"])

#makes chart that summarizes the data
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)

hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")

legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max)) 
legend("topright", legend_text)

#extract gene names and transcript names
gene_names=data.frame(SYMBOL=unique(rownames(gene_bg_chrx)))

## create sample meta data frame
phenotype_table= data.frame(id=sampleNames(bg_chrX), group=rep(c("Female","Male"), each=3))

pData(bg_chrX) =pheno_data

#plot transcripts for 1 gene for 1 sample
plotTranscripts(gene='MSTRG.6938', bg_chrX, samples='male1_lane7', meas='FPKM', colorby='transcript', main='transcripts from gene MSTRG.6938: male1_lane7, FPKM')

plotMeans('MSTRG.6938',bg_chrX,groupvar='sex' ,meas='FPKM', colorby='transcript')

## differential transcript expression
results_txns = stattest(bg_chrX, feature='transcript', getFC = T, covariate='sex',meas='FPKM' )

#add names and IDs
t.ids=Whole_bg_chrx[,c(1,6)]
t_names=unique(Whole_bg_chrx[,c(1,6)])

results_txns_merged = merge(results_txns,t.ids,by.x=c("id"),by.y=c("t_id"))
head(results_txns_merged)

# Calculate differentially expressed genes and use FPKM in calculating # # # # differential gene expression
results_genes = stattest(bg_chrX, feature="gene", covariate="sex", getFC=TRUE, meas="FPKM")

## Compare the data before and after normalization. boxplot with and without log transformation
par(mfrow=c(1,2))

boxplot(gene_expression, col=rainbow(6),  las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 6 samples")

boxplot(log2(gene_expression+1), col= Beach,  las=2, ylab="log2(FPKM)", main="log transformed distribution of FPKMs for all 6 samples")
dev.off()

## FPKM values are not logged. Hence fold change (FC) is not also logged. Log # fold changes and store it in logfc columnresults_genes[,"logfc"] = log2(results_genes[,"fc"])
results_genes[,"logfc"] = log2(results_genes[,"fc"])

# Identify the genes (rows) with adjusted p-value (i.eq-value) < 0.05 
qsig=which(results_genes$qval<0.05)

# draw histogram
hist(results_genes[qsig,"logfc"], breaks=50, col="seagreen", xlab="log2(Fold change) male vs female", main="Distribution of differential expression values")
abline(v=c(-2,2), col="black", lwd=2, lty=2)
legend("topright", "Fold change >2 and <-2", lwd=2, lty=2)

# Convert the matrix to data
gene_expression=as.data.frame(gene_expression)
write.table(gene_expression, "gene_expression.txt", sep="\t")

gene_expression$female=rowMeans(gene_expression[, c(1:8)])
gene_expression$male=rowMeans(gene_expression[,c(9:12)])

#to avoid 0
x=log2(gene_expression[,"female"]+1)
y=log2(gene_expression[,"male"]+1)

plot(x=x, y=y, pch=1, cex=2, xlab="Female FPKM (log2)", ylab="Male (log2)", main="Male vs Female FPKMs")
abline(a=0, b=1)

xqsig=x[qsig]
yqsig=y[qsig]
points(x=xqsig, y=yqsig, col="green", pch=19, cex=2)

fsig=which(abs(results_genes$logfc)>4)
xfsig=x[fsig]
yfsig=y[fsig]
points(x=xfsig, y=yfsig, col="red", pch=1, cex=2)

legend_text = c("Significant by Q value", "Significant by Fold change")
legend("topright", legend_text,bty="n",pch = c(19,19), col=c("green","red"))

# label the significant genes
textxy(xfsig,yfsig, cex=0.8, labs=row.names(gene_expression[fsig,]))
# add red line through 0
abline(v=0, col="red", lwd=3)

# add red line through fold change 4 (log2,2)
abline(v=c(4,-4), col="red", lwd=3)
abline(h=c(-4,4), col="red",lwd=3)

#volcano plot
# Filter genes by log fold change by 16 fold
fc_sig_results_genes=which(abs(results_genes$logfc)>4)

fc_sig_results_genes_plot=results_genes[fc_sig_results_genes,]  

plot(results_genes$logfc,results_genes$qval, col="steelblue", pch=1) 

#abline
abline(v=c(2,-2), col="red", lwd=3)
abline(h=0.05, col="red",lwd=3)

# highlight the genes with color  
points(fc_sig_results_genes_plot$logfc,fc_sig_results_genes_plot$qval, col="green", pch=16) 

# label the significant genes
textxy(fc_sig_results_genes_plot$logfc,fc_sig_results_genes_plot$qval, labs=fc_sig_results_genes_plot$id, cex=1.2)

colors = colorRampPalette(c("white", "blue","red","green","yellow"))
par(mfrow=c(1,2))
plot(x,y)
smoothScatter(x,y, colramp = colors)

# Identify the genes (rows) below p-value 0.05
sigpi = which(results_genes[,"pval"]<0.05)

# Extract p-significant genes in a separate object
sigp = results_genes[sigpi,]

sigde = which(abs(sigp[,"logfc"]) >= 2)

# Extract and store the statistically significant genes (rows) that are upregulated/ downregulated by 4 fold
sig_tn_de = sigp[sigde,]

# Order by q value, followed by differential expression
sorted_sig_tn_de = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"logfc"]), decreasing=FALSE)

output = sig_tn_de[sorted_sig_tn_de,c("id","fc","pval","qval","logfc")]
write.table(output, file="ballgown/SigDE.txt", sep="\t", row.names=FALSE, quote=FALSE)

sig_gene_expression=gene_expression[rownames(gene_expression) %in% sig_tn_de$id,]

#remove female and male columns
sig_gene_expression=sig_gene_expression[,-c(7:9)]

# for pheatmap function, column names and row names of data and pdata mush be identical# change the row names
rownames(pheno_data)=pheno_data[,1]

# remove the id column
phenotype_table=subset(pheno_data, select = -c(id) )

# change the colnames to match with the sample names
colnames(sig_gene_expression)=row.names(pheno_data)

library(pheatmap)
pheatmap(as.matrix(sig_gene_expression), scale = "row", clustering_distance_rows = "correlation", clustering_method = "complete",annotation_col = pheno_data , main="Significant genes",fontsize_col=14, fontsize_row = 6 ,color = c("green","red"))



#PcA plot
pca_data=prcomp(t(sig_gene_expression))

# Calculate PCA component percentages
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)

df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = colnames(sig_gene_expression), condition = rep(c("Male","Female")))

library(ggplot2)
library(ggrepel)

ggplot(df_pca_data, aes(PC1,PC2, color = sample))+
  geom_point(size=8)+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")")) 

ggplot(df_pca_data, aes(PC1,PC2, color = condition))+
  geom_point(size=8)+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))+
  geom_text_repel(aes(label=sample),point.padding = 0.75)
