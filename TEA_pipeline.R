### module load gcc/6.2.0
### module load R/3.3.0


# .libPaths("/mnt/work1/users/bhklab/Rlib")
# source("TEA_pipeline.R")

########################
### clean workspace
rm(list = ls())

########################
### load libraries
require(gdata)
require(e1071)
require(genefu)
library(Biobase)
require(xtable)
library(biomaRt)
library(gplots)
library(devtools)
library(preprocessCore)
library(rgl)
library(qpcR)
library(data.table)
library(piano)
library(snowfall)
library(utils)
library(Hmisc)
library(WriteXLS)
# devtools::install_github(repo="bhklab/PharmacoGx")
library(PharmacoGx)
library(survcomp)

########################
### global parameters

### set global R options
options(stringsAsFactors=FALSE)

### create directory for all the results and intermediate files
OutDir <- file.path("Output")
dir.create(OutDir, showWarnings=FALSE, recursive=TRUE)

### create directory for enrichment results
GSEADir <- file.path(OutDir, "Drug_Tissue_Associations")
dir.create(GSEADir, showWarnings=FALSE, recursive=TRUE)

### quantile for AUC in each tissue type of interest
quantileAUC <- 0.75

### Min and max number of cell lines in a tissue type
TissueSize <- c(15, 200)

FDRcutoff <- 0.05

### number of permutations for enrichment analysis
nperm <- 1000000

### number of CPU cores used for parallelization, use NULL for all the cores minus one
nbcore <- NULL
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore - 1 }
options("mc.cores"=nbcore)

########################
### load list of in vitro drug screening datasets
PsetList <- list("CCLE", "gCSI", "CTRPv2", "GDSC1000")
myfn <- file.path(OutDir, "psets.rds")
if (!file.exists(myfn)) {
    PsetVec <- lapply(PsetList, PharmacoGx::downloadPSet, saveDir=OutDir)
    names(PsetVec) <- PsetList
    saveRDS(object=PsetVec, file=myfn)
} else {
    PsetVec <- readRDS(file=myfn)
}

########################
### update the PSets for tissue enrichment analysis

### split lung into NSCLC and SCLC
NSCLC_cellines <- PsetVec$GDSC1000@cell$Sample.Name[which(PsetVec$GDSC1000@cell$GDSC.Tissue.descriptor.1 == "lung_NSCLC" & PsetVec$GDSC1000@cell$tissueid == "lung")]
SCLC_cellines <- PsetVec$GDSC1000@cell$Sample.Name[which(PsetVec$GDSC1000@cell$GDSC.Tissue.descriptor.1 == "lung_SCLC" & PsetVec$GDSC1000@cell$tissueid == "lung")]

### CCLE
PsetVec$CCLE@cell <- cbind(PsetVec$CCLE@cell, "tissueid_TEA"=PsetVec$CCLE@cell[ , "tissueid"])
PsetVec$CCLE@cell[ , "tissueid_TEA"] <- as.character(PsetVec$CCLE@cell[ , "tissueid_TEA"])
PsetVec$CCLE@cell[!is.na(PsetVec$CCLE@cell[ , "tissueid_TEA"]) & PsetVec$CCLE@cell[ , "tissueid_TEA"] == "", "tissueid_TEA"] <- NA
PsetVec$CCLE@cell[!is.na(PsetVec$CCLE@cell[ , "tissueid_TEA"]) & PsetVec$CCLE@cell[ , "tissueid_TEA"] == "lung", "tissueid_TEA"] <- NA
PsetVec$CCLE@cell[which((PsetVec$CCLE@cell[ ,"Hist.Subtype1"] == "adenocarcinoma" |
    PsetVec$CCLE@cell[ ,"Hist.Subtype1"] == "non_small_cell_carcinoma" |
    PsetVec$CCLE@cell[ ,"Hist.Subtype1"] == "squamous_cell_carcinoma") & 
    PsetVec$CCLE@cell[ ,"tissueid"] == "lung"), "tissueid_TEA"] <- "NSCLC"
PsetVec$CCLE@cell[which(PsetVec$CCLE@cell[ ,"Hist.Subtype1"] == "small_cell_carcinoma" & 
    PsetVec$CCLE@cell[ ,"tissueid"] == "lung"),"tissueid_TEA"] <- "SCLC"

### gCSI
PsetVec$gCSI@cell <- cbind(PsetVec$gCSI@cell, "tissueid_TEA"=PsetVec$gCSI@cell[ , "tissueid"])
PsetVec$gCSI@cell[ , "tissueid_TEA"] <- as.character(PsetVec$gCSI@cell[ , "tissueid_TEA"])
PsetVec$gCSI@cell[!is.na(PsetVec$gCSI@cell[ , "tissueid_TEA"]) & PsetVec$gCSI@cell[ , "tissueid_TEA"] == "", "tissueid_TEA"] <- NA
PsetVec$gCSI@cell[!is.na(PsetVec$gCSI@cell[ , "tissueid_TEA"]) & PsetVec$gCSI@cell[ , "tissueid_TEA"] == "lung", "tissueid_TEA"] <- NA
PsetVec$gCSI@cell[which(PsetVec$gCSI@cell[ ,"CellLineName"] %in% NSCLC_cellines & PsetVec$gCSI@cell[ ,"tissueid"] == "lung"), "tissueid_TEA"] <- "NSCLC"
PsetVec$gCSI@cell[which(PsetVec$gCSI@cell[ ,"CellLineName"] %in% SCLC_cellines & PsetVec$gCSI@cell[,"tissueid"] == "lung"), "tissueid_TEA"] <- "SCLC"

### CTRPv2
PsetVec$CTRPv2@cell <- cbind(PsetVec$CTRPv2@cell, "tissueid_TEA"=PsetVec$CTRPv2@cell[ , "tissueid"])
PsetVec$CTRPv2@cell[ , "tissueid_TEA"] <- as.character(PsetVec$CTRPv2@cell[ , "tissueid_TEA"])
PsetVec$CTRPv2@cell[!is.na(PsetVec$CTRPv2@cell[ , "tissueid_TEA"]) & PsetVec$CTRPv2@cell[ , "tissueid_TEA"] == "", "tissueid_TEA"] <- NA
PsetVec$CTRPv2@cell[!is.na(PsetVec$CTRPv2@cell[ , "tissueid_TEA"]) & PsetVec$CTRPv2@cell[ , "tissueid_TEA"] == "lung", "tissueid_TEA"] <- NA
PsetVec$CTRPv2@cell[which((PsetVec$CTRPv2@cell[,"ccle_hist_subtype_1"] == "adenocarcinoma" |
    PsetVec$CTRPv2@cell[ , "ccle_hist_subtype_1"] == "non_small_cell_carcinoma" |
    PsetVec$CTRPv2@cell[ , "ccle_hist_subtype_1"] == "squamous_cell_carcinoma") & 
    PsetVec$CTRPv2@cell[ , "tissueid"] == "lung"), "tissueid_TEA"] <- "NSCLC"
PsetVec$CTRPv2@cell[which(PsetVec$CTRPv2@cell[ , "ccle_hist_subtype_1"] == "small_cell_carcinoma" & 
    PsetVec$CTRPv2@cell[ ,"tissueid"] == "lung"), "tissueid_TEA"] <- "SCLC"

### GDSC1000
PsetVec$GDSC1000@cell <- cbind(PsetVec$GDSC1000@cell, "tissueid_TEA"=PsetVec$GDSC1000@cell[ , "tissueid"])
PsetVec$GDSC1000@cell[ , "tissueid_TEA"] <- as.character(PsetVec$GDSC1000@cell[ , "tissueid_TEA"])
PsetVec$GDSC1000@cell[!is.na(PsetVec$GDSC1000@cell[ , "tissueid_TEA"]) & PsetVec$GDSC1000@cell[ , "tissueid_TEA"] == "", "tissueid_TEA"] <- NA
PsetVec$GDSC1000@cell[!is.na(PsetVec$GDSC1000@cell[ , "tissueid_TEA"]) & PsetVec$GDSC1000@cell[ , "tissueid_TEA"] == "lung", "tissueid_TEA"] <- NA
PsetVec$GDSC1000@cell[which(PsetVec$GDSC1000@cell[,"GDSC.Tissue.descriptor.1"] == "lung_NSCLC" & PsetVec$GDSC1000@cell[ ,"tissueid"] == "lung"), "tissueid_TEA"] <- "NSCLC"
PsetVec$GDSC1000@cell[which(PsetVec$GDSC1000@cell[,"GDSC.Tissue.descriptor.1"] == "lung_SCLC" & PsetVec$GDSC1000@cell[ ,"tissueid"] == "lung"), "tissueid_TEA"] <- "SCLC"

########################
### list drugs and cell lines and drugs in all datasets
### drugs
drugs <- sapply(PsetVec, PharmacoGx::drugNames)
drugs <- sort(unique(do.call(c, drugs)))
drugsMat <- matrix(NA, nrow=length(drugs), ncol=length(PsetVec), dimnames=list(drugs, names(PsetVec)))
for (PsetIter in 1:length(PsetVec)) {
  drugsMat[PharmacoGx::drugNames(PsetVec[[PsetIter]]), names(PsetVec)[PsetIter]] <- "YES"
}
### cell lines
cells <- sapply(PsetVec, PharmacoGx::cellNames)
cells <- sort(unique(do.call(c, cells)))
cellsMat <- matrix(NA, nrow=length(cells), ncol=length(PsetVec), dimnames=list(cells, names(PsetVec)))
for (PsetIter in 1:length(PsetVec)) {
  cellsMat[PharmacoGx::cellNames(PsetVec[[PsetIter]]), names(PsetVec)[PsetIter]] <- "YES"
}
### tissues
tissues <- sapply(PsetVec, function(x) { return (sort(unique(PharmacoGx::cellInfo(x)[ , "tissueid_TEA"]))) })
tissues <- sort(unique(do.call(c, tissues)))
tissuesMat <- matrix(NA, nrow=length(tissues), ncol=length(PsetVec), dimnames=list(tissues, names(PsetVec)))
for (PsetIter in 1:length(PsetVec)) {
  tt <- sort(unique(PharmacoGx::cellInfo(PsetVec[[PsetIter]])[ , "tissueid_TEA"]))
  tissuesMat[tt, names(PsetVec)[PsetIter]] <- "YES"
}
### save
ll <- list("Drugs"=data.frame(drugsMat), "Cell lines"=data.frame(cellsMat), "Tissue Types"=data.frame(tissuesMat))

for(ll_c in 1:length(ll)){
  xlsx::write.xlsx(ll[[ll_c]], file = file.path(OutDir, sprintf("Dataset_Info.xlsx")), row.names = TRUE, append = TRUE, sheetName = names(ll)[ll_c], showNA = FALSE)
}

#WriteXLS::WriteXLS("ll", ExcelFileName=file.path(OutDir, sprintf("Dataset_Info.xlsx")), row.names=TRUE)

########################

### run tissue enrichment analysis with original AUC
Adjustment <- FALSE
source("TEA_analysis.R")

### run tissue enrichment analysis with AUC values adjusted for genel level of drug sensitivity
Adjustment <- TRUE
source("TEA_analysis.R")

### combine the results
source("TEA_postprocess.R")

### compute predictability of significant interactions
source("TEA_predictability.R")

source("TEA_figure_generator.R")

########################
### save session info

write(toLatex(sessionInfo(), locale=FALSE), file=file.path(OutDir, "sessionInfoR.tex"), append=FALSE)

### end

