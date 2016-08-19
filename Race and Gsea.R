#Packages Required#

must.install <- FALSE
if(require(RTCGAToolbox)){
  si <- devtools::session_info()
  if(!grepl("link-ny", si$packages[si$packages[, 1] == "RTCGAToolbox", 5], ignore.case = TRUE)){
    must.install <- TRUE
  }
  if(!require("TCGAmisc")){
    must.install <- TRUE
  }
}else{
  must.install <- TRUE
}
if(must.install){
  BiocInstaller::biocLite(c("Link-NY/RTCGAToolbox", "waldronlab/TCGAmisc"))
}


library("TCGAmisc")
library("gage")
library("readr")
library(RTCGAToolbox)

if(file.exists("prad_eset.rds" & file.exists("prad_mirna_eset.rds"))){
  prad_eset <- readRDS("prad_eset.rds")
  prad_mirna_eset <- readRDS("prad_mirna_eset.rds")
}else{
  (rundates <- getFirehoseRunningDates())
  prad <- getFirehoseData("PRAD", runDate="20151101",
                        miRNASeq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE)
  prad_eset <- extract(prad, "rnaseq2genenorm")
  prad_mirna_eset <- extract(prad, "miRNASeq_Gene")
  saveRDS(prad_eset, file="prad_eset.rds")
  saveRDS(prad_mirna_eset, file="prad_mirna_eset.rds")
}


# 1304 targets and Race Variable upload#



miRNA1304.targets <- read_csv("https://raw.githubusercontent.com/lwaldron/tcga_prad/master/miRNA1304%20targets.csv")
racevar <- read_csv("https://raw.githubusercontent.com/lwaldron/tcga_prad/master/racevariable.csv")


#Matching Race Variable from file to the expr(eset) race variable and creating two separate classes (column vectors) for blacks and whites#


racevar$patientID <- barcode(racevar$X)
esetID <- barcode(colnames(exprs(prad_eset)))

racevar <- racevar[match(esetID, racevar$patientID), ]

if( identical(esetID, racevar$patientID) ){
  prad_eset$race <- racevar$race
}else{
  prad_eset$race <- NA
}

blacks <- which(prad_eset$race == "black or african american")
whites <- which(prad_eset$race == "white")

## kegg.db uses Entrez Gene identifiers and prad_eset uses gene symbols, 
## so these need to be mapped.  An easier solution would be to get the mSigDB 
## gene sets that already use gene symbols.

exdat <- log2(exprs(prad_eset) + 1)
library(org.Hs.eg.db)
select(org.Hs.eg.db, head(rownames(exdat)), "ENTREZID")
symbols <- mapIds(org.Hs.eg.db, rownames(exdat), column="ENTREZID", keytype="SYMBOL")
mapped <- !is.na(symbols)
exdat <- exdat[mapped, ]
rownames(exdat) <- symbols[mapped]

#Gage Analysis for Race Response#

kegg.p <- gage(exprs = exdat, gsets = kegg.gs, ref= whites, samp = blacks, compare="unpaired")
kegg.sig <- sigGeneSet(gse16873.kegg.p, outname="kegg")

