library(gage)
library(RTCGAToolbox)
rundates <- getFirehoseRunningDates()
analysisdates <- getFirehoseAnalyzeDates()
prad <- getFirehoseData("PRAD", runDate=rundates[1],
                        gistic2_Date=analysisdates[1], RNAseq_Gene=TRUE, 
                        miRNASeq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE, CNA_SNP=FALSE,
                        CNV_SNP=FALSE, CNA_Seq=FALSE, CNA_CGH=FALSE,  Methylation=FALSE,
                        Mutation=FALSE, mRNA_Array=TRUE, miRNA_Array=TRUE, RPPA=FALSE)

eset <- extract(prad, "rnaseq2genenorm")
eset1<- extract(prad, "miRNASeq_Gene")

y<-log2(exprs(eset)+1)

zz <- exprs(eset1)["hsa-mir-1304",]
z <- log2(zz+1)


#Sanity Check, data preparation#
gagePrep(exprs= y, samp = z, same.dir = TRUE, rank.test = FALSE, use.fold = TRUE, weights = NULL)

targets <- read.table("miRNA1304 targets.csv", quote="\"", comment.char="", as.is=TRUE)
colnames(targets) <- "miRNA1304 targets"

#Gage Results(miRNA-1304 response) for greater than, less than, and stat geometric mean#
results<-gage(exprs = y, gsets = as.list(targets), samp = z, full.table = TRUE, saaPrep = gagePrep, saaTest = gs.tTest, same.dir = TRUE, test4up=TRUE)


#Essenial Genes from Target List#

essGene(targets, exprs= y, ref = NULL, samp = z, use.fold = TRUE, rank.abs = FALSE, use.chi = FALSE, chi.p = 0.05)

#Attempt to get essential gene list although I do not understand how to get the p-value matrix from results to setp#
esset.grp(setp = , exprs=y, gsets = list(targets) , ref = NULL, samp = z, test4up = TRUE, same.dir = TRUE, use.fold = TRUE, cutoff = 0.05, use.q = FALSE, pc = 10^-10, output = TRUE, outname = "esset.grp", make.plot = TRUE, pdf.size = c(7, 7), core.counts = FALSE, get.essets = TRUE, bins = 10, bsize = 1, cex = 0.5, layoutType = "circo", name.str =c(10, 100))

#Attempt to make scatterplot or heatmap of essential genes from targets, although I could not get this to work correctly due to length issues with exprs and samp#
geneData(targets1, exprs = y, samp = z, outname = "array", txt = TRUE, heatmap = FALSE, scatterplot = TRUE, samp.mean = FALSE,pdf.size = c(7, 7), cols = NULL, scale = "row", limit = NULL, label.groups = TRUE)