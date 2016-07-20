## RNA-seq analysis with DESeq2
source('http://bioconductor.org/biocLite.R')
#biocLite('DESeq2')
library('DESeq2')


# Import & pre-process ----------------------------------------------------

# Import data from HTseqcount Call
directory<- # Directory of HTseq count output (must be .txt)
#use grep to search for the 'htseq' part of filename to collect files
sampleFiles<-grep('htseq',list.files(directory),value=TRUE)
# view sampleFiles
sampleFiles

sampleCondition<-c('N2_WT','eat2_mut','nhr62_mut','e2_n6')
sampleBatch <- c("Batch1","Batch1","Batch1","Batch1","Batch2","Batch2","Batch2","Batch2","Batch3","Batch3","Batch3","Batch3")

sampleCondition <- factor(sampleCondition, levels=c('N2_WT','eat2_mut','nhr62_mut','e2_n6'))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, genotype=sampleCondition ,Batch = sampleBatch)
sampleTable 

# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)

## view ddsHTseq - should give summary of class, data, etc.
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=formula(~ Batch + genotype))
ddsHTSeq

#Running the differential expression pipeline
dds<-DESeq(ddsHTSeq)
plotDispEsts(dds)
show(resultsNames(dds))



# Get differential expression results
res1=results(dds,contrast=c("genotype","N2_WT","eat2_mut"), alpha=0.05)

# order results by log2FoldChange (most significant to least) and extract significant 
res1 <- res1[order(res1$log2FoldChange),]
WT_vs_eat2_mut <- res1[ which(res1$padj < 0.05), ] 

#write the table to a csv file
write.csv(as.data.frame(WT_vs_eat2_mut),file='DE_gene_WT_vs_eat2_mut.csv')

# MA plot of RNAseq data for entire dataset
plotMA(res1,main='DESeq2 N2_WT vs eat2_mut',ylim=c(-3,3))
dev.copy(png,'deseq2_MAplot_res1.png')
dev.off()

# Make a basic volcano plot
with(res1, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot N2_WT vs eat2_mut", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res1, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res1, abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res1, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green"))


#### -------------- Regularized log transformation for clustering/heatmaps, etc --------------- ####
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

#clustering analysis
library('RColorBrewer')
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]

# heatmap of data
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:1000]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(assay(rld)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="TITLE")
dev.copy(png, "DATE-DESeq2-HEATMAP.png")
dev.off()

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(genotype,sampleFiles , sep=' : '))
hc <- hclust(distsRL)
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmaptest<- heatmap.2(mat, Rowv=as.dendrogram(hc),
                        symm=TRUE, trace='none',
                        col = rev(hmcol), margin=c(20, 20))
dev.copy(png,'deseq2_heatmaps_samplebysample.png')
dev.off()

#Principal components plot
print(plotPCA(rld, intgroup=c('genotype')))
print(plotPCA(rld,intgroup=c("genotype", "Batch")))
dev.copy(png,'deseq2_pca.png')
dev.off()


###          GO Enrichment ANlaysis   ###


library("goseq")
library("geneLenDataBase")
library("org.Ce.eg.db")
library("GenomicFeatures")

#https://wikis.utexas.edu/display/bioiteam/GO+Enrichment+using+goseq

ALL_N2vseat2<-rownames(res1)
ALL_N2vseat2_vector<-t(as.vector(ALL_N2vseat2))


WT_vs_eat2_mut<-rownames(WT_vs_eat2_mut)
WT_vs_eat2_mut<-t(as.vector(WT_vs_eat2_mut))

#F) Construct a new vector that adds a 0 next to every gene that is not in our DEG list and a 1 next to every gene that is in our DEG list.

gene_vector_N2vseat2=as.integer(ALL_N2vseat2_vector%in%WT_vs_eat2_mut)
names(gene_vector_N2vseat2)=ALL_N2vseat2_vector

#Weigh the gene vector by lengths of our genes. This length information is pulled from a different package and we need to specifying the exact name of our genome and gene IDs in that package. 
pwf2=nullp(gene_vector_N2vseat2,"ce11","ensGene")

#If the information on  median transcript length is not provided you have to create it (Necessary in our case)
#txdb = makeTxDbFromBiomart(biomart = 'ensembl', dataset = 'celegans_gene_ensembl')
#txdbByGene = transcriptsBy(txdb, 'gene')
#lengthData = median(width(txdbByGene))
#pwf2=nullp(gene_vector_N2vseat2,"ce11","ensGene",bias.data=lengthData)



#Find enriched GO terms
GO.wall=goseq(pwf2,"ce11","ensGene",use_genes_without_cat=TRUE)

#How many enriched GO terms do we have
class(GO.wall)
head(GO.wall)
nrow(GO.wall)

#If you wish to identify categories significantly enriched/unenriched below some p-value cutoff, it is necessary to first apply some kind of multiple hypothesis testing correction.
GO.wall$padj=p.adjust(GO.wall$over_represented_pvalue,method="BH")
enriched.GO<-GO.wall[ which(GO.wall$padj<0.05),]
#write the table to a csv file
write.csv(as.data.frame(enriched.GO),file='Go_enrichment.csv',row.names = F)

#PLot Gene Ontology Tree with RamiGO
library("RamiGO")
ramigo.A = getAmigoTree (goIDs = enriched.GO$category, pvalues = enriched.GO$padj, pcolors = c('white', 'magenta'), psplit = c(0.1, 0.05, 0.025, 0.01, 0.00025), filename = 'RamiGO_A', picType = 'png', modeType = 'amigo', saveResult = TRUE)

