######  Converting Refseq identifiers to HGNC symbols with biomaRt ###########

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
library(readr)


racevar <- read_csv("https://raw.githubusercontent.com/pchernek/tcga_prad/master/miRDB_v5.0_prediction_result.csv")

x<-miRDB_v5.0_prediction_result

x<-as.vector(x)

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

listFilters(mart)

ann <- getBM(c("hgnc_symbol"),"refseq_mrna_predicted", x[1], mart)

ann1 <- getBM(c("hgnc_symbol"),"refseq_mrna", x[1], mart)
