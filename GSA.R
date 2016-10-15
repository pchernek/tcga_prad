library(readr)

library(dplyr)
library(BiocInterfaces)
library(RTCGAToolbox)
library(GSA)

rundates <- getFirehoseRunningDates()



analysisdates <- getFirehoseAnalyzeDates()
prad <- getFirehoseData("PRAD", runDate=rundates[1],
                        gistic2_Date=analysisdates[1], RNAseq_Gene=FALSE, 
                        miRNASeq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE, CNA_SNP=FALSE,
                        CNV_SNP=FALSE, CNA_Seq=FALSE, CNA_CGH=FALSE,  Methylation=FALSE,
                        Mutation=FALSE, mRNA_Array=FALSE, miRNA_Array=FALSE, RPPA=FALSE)

choices <- tolower(gsub("_", "", c("RNAseq_Gene", "miRNASeq_Gene",
                                   "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
                                   "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
                                   "miRNA_Array", "RPPA")))

dses <- lapply(choices, function(choice) try(extract(prad, choice, 
                                                     clinic=TRUE),
                                             silent=TRUE))
names(dses) <- choices
dses


eset<- TCGAextract(prad, "RNAseq2_Gene_Norm")



colnames(eset)<-TCGAbarcode(colnames(eset), sample = TRUE)

eset<-as.matrix(eset)
tumors<-grepl("01", colnames(eset))





eset<-eset[, tumors]


z3= substr(colnames(eset),1,nchar(colnames(eset))-4)

z3<-tolower(colnames(eset))


z3= substr(colnames(eset),1,nchar(colnames(eset))-4)


colnames(eset)<-z3


racevar <- read_csv("https://raw.githubusercontent.com/lwaldron/tcga_prad/master/racevariable.csv")

racevar$X <-(gsub("\\.", "-", racevar$X))
racevar$X
racevar$patientID <- TCGAbarcode(racevar$X)



eset<-eset[, match(racevar$patientID, z3)]



Race <- z3[match(z3, racevar$X)]


Race
Race<-as.data.frame(Race)
names(Race)[1]<-"V2"


newRace<- left_join(Race, racevar, by = c("V2" = "X"))


newRace<-as.data.frame(newRace)


Clinical1<-prad@Clinical[, c(2, 8, 9, 16, 17, 18)]




x<-setDT(Clinical1, keep.rownames = TRUE)
x<-as.data.frame(Clinical1)



x$rn<-(gsub("\\.", "-", (x$rn)))




masterRace <- newRace$V2[match(x$rn,newRace$V2)]
masterRace<-masterRace[!duplicated(masterRace)]

masterRace<-as.data.frame(masterRace)
names(masterRace)[1]<-"V3"


p<-left_join(masterRace, newRace, by = c("V3" = "V2"))
y<-left_join(p, x, by = c("V3" = "rn") )


eset<-eset[, match(y$V3, z3)]

y= subset(y, y[, 2] %in% c("white", "black or african american"))
esetFINAL<-eset[, match(y$V3, colnames(eset))]








##GSA Analysis RNA seq and Race Response##

x2<-which((y[, 2])=='white')


x3<-which(y[, 2]=='black or african american')


subset1<-log2((esetFINAL)+ 1)

subset1<-as.matrix(subset1)




x<-subset1



            
          
w<-c(x3, x2)


genenames<- rownames(x)
genenames<-as.vector(genenames)
genenames



geneset.obj<- GSA.read.gmt("https://raw.githubusercontent.com/pchernek/tcga_prad/master/msigdb.v5.2.symbols.gmt.txt")


GSA.obj<-GSA(x,w, genenames = genenames,  genesets=geneset.obj$genesets, nperms=100)


GSA.listsets(GSA.obj, FDRcut = 0.05, geneset.names = geneset.obj$geneset.names)

GSA.plot(GSA.obj, fac = 0.05,  FDRcut = .95)




#