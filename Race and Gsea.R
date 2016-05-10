#Packages Required#



library("TCGAmisc")
library("gage")
library("readr")


# Data types needed from Prad#

eset <- extract(prad, "rnaseq2genenorm")

eset1<- extract(prad, "miRNASeq_Gene")


# mRNA and miRNA 1304 Expression Datasets#

y<-log2(exprs(eset)+1)

zz <- exprs(eset1)["hsa-mir-1304",]

z <- log2(zz+1)


# 1304 targets and Race Variable upload#



miRNA1304.targets <- read_csv("https://raw.githubusercontent.com/lwaldron/tcga_prad/master/miRNA1304%20targets.csv")
racevar <- read_csv("https://raw.githubusercontent.com/lwaldron/tcga_prad/master/racevariable.csv")


#Matching Race Variable from file to the expr(eset) race variable and creating two separate classes (column vectors) for blacks and whites#


racevar$patientID <- barcode(racevar$X)
esetID <- barcode(colnames(exprs(eset)))
esetRace <- racevar$race[match(esetID, racevar$patientID)]
blacks <- which(esetRace == "black or african american")
whites <- which(esetRace == "white")

blacks

whites



#Gage Analysis for Race Response#

gagePrep(exprs= y, samp = blacks, same.dir = TRUE, rank.test = FALSE, use.fold = TRUE, weights = NULL)

gage(exprs = y, gsets = as.list(miRNA1304.targets), ref= whites , samp = blacks, full.table = TRUE, saaPrep = gagePrep, saaTest = gs.tTest, same.dir = TRUE, test4up=TRUE, compare = "unpaired")