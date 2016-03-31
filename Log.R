source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("liNk-NY/RTCGAToolbox")
biocLite("devtools")
library(RTCGAToolbox)
rundates <- getFirehoseRunningDates()
analysisdates <- getFirehoseAnalyzeDates()
prad <- getFirehoseData("PRAD", runDate=rundates[1],
                        gistic2_Date=analysisdates[1], RNAseq_Gene=TRUE, 
                        miRNASeq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE, CNA_SNP=TRUE,
                        CNV_SNP=TRUE, CNA_Seq=TRUE, CNA_CGH=TRUE,  Methylation=TRUE,
                        Mutation=TRUE, mRNA_Array=TRUE, miRNA_Array=TRUE, RPPA=TRUE)

choices <- tolower(gsub("_", "", c("RNAseq_Gene", "miRNASeq_Gene",
                                   "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
                                   "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
                                   "miRNA_Array", "RPPA")))

dses <- lapply(choices, function(choice) try(extract(prad, choice, 
                                                     clinic=TRUE),
                                             silent=TRUE))
names(dses) <- choices
dses

eset <- extract(prad, "rnaseq2genenorm")
eset1<- extract(prad, "miRNASeq_Gene")

x<-exprs(eset)
y<-log(x+1)
write.csv(y, file = "log+1 of messengerRNA Expression")
w<-exprs(eset1)[156,]
View(w)
View(x)
z<-log(w)
write.csv(z, file = "log of micro1304 expression")