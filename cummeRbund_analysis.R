library(cummeRbund)
library(cairoDevice)
setwd()<-# path/to/cuffdiff_output/
  cuff_data <-readCufflinks(rebuild = TRUE)
cuff_data 

#Quality Assessment of data-Evaluating model fit
csDensity(genes(cuff_data)) #Plot the distribution of expression levels for each sample

dispersionPlot(genes(cuff_data))

# Identifying outlier replicates
pBoxRep<-csBoxplot(genes(cuff_data),replicates=T)
pBoxRep

pDendro<-csDendro(genes(cuff_data),replicates=T)
pDendro

#Boxplot View of Expression Levels - distributions of FPKM 

genes.scv<-fpkmSCVPlot(genes(cuff_data))
genes.scv
isoforms.scv<-fpkmSCVPlot(isoforms(cuff_data))
isoforms.scv


#A matrix of pairwise scatterplots can be drawn using the csScatterMatrix() method.
s<-csScatterMatrix(genes(cuff_data))
s

#Compare the expression of each gene in two conditions with a scatter plot
scatterN2_WT_eat2_mut <- csScatter(genes(cuff_data), 'N2_WT','eat2_mut') 
scatterN2_WT_eat2_mut 

#Create a volcano plot to inspect differentially expressed genes
volcanoN2_WT_eat2_mut <- csVolcano(genes(cuff_data),'N2_WT','eat2_mut',alpha=0.05, showSignificant=T)#ShowSignificant to show in red genes DE
volcanoN2_WT_eat2_mut

####### Record differentially expressed genes to files for use in downstream analysis ########

#	retrieve	the	gene-level	differential	expression	data	
gene_diff_data	=	diffData(genes(cuff_data))
#	how	many	‘genes’ are	there?
nrow(gene_diff_data)
#	from	the	gene-level	differential	expression	data,	extract	those	that are	labeled	as	significantly	different.
sig_gene_data	=	subset(gene_diff_data,(significant=='yes'),row.names=T)
#	how	many genes	are	significantly	DE	according	to	these	criteria?
nrow(sig_gene_data)
#	Examine	the	entries	at	the	top	of	the	unsorted	data	table:
head(sig_gene_data)
# Record differentially expressed genes by group
N2_vs_eat2<-subset(sig_gene_data, sample_1=='N2_WT' & sample_2=='eat2_mut')
#	Write	the	list	of	significantly	differentially	expressed	genes	to	a	file	like	so:
write.table(N2_vs_eat2,file='sig_diff_genes.txt',sep	='\t',	quote	=	F)

####               Inspect the Differentially Expressed Transcripts                 ####
isoform_diff_data_WT_eat2_mut <- diffData(isoforms(cuff_data),'N2_WT','eat2_mut')
sig_isoform_data_WT_eat2_mut <- subset(isoform_diff_data_WT_eat2_mut, (significant == 'yes'))
nrow(sig_isoform_data_WT_eat2_mut)

write.table(sig_isoform_data_WT_eat2_mut,file='sig_isoform_data_WT_eat2_mut.txt',sep	='\t',	quote	=	F)


#Plot expression levels for genes of interest with bar plots 
mygene <- getGene(cuff_data, 'XLOC_000040') #	get that	gene	‘object’	from	cummeRbund
mygene
expressionBarplot(mygene,logMode=T) 	#	plot	the	expression	values	for	the	gene	under	each	condition
expressionBarplot(isoforms (mygene),replicates=T) #Plot individual isoform expression levels of selected genes of interest with bar plots :

head(fpkm(isoforms(mygene)))
gl<-expressionPlot(mygene)
gl
