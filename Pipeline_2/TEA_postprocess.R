#################
rm(list = ls())
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
# devtools::install_github(repo="bhklab/PharmacoGx")
library(PharmacoGx)
###################
Adjustment <- TRUE
##############################
CTRPv2 <- downloadPSet("CTRPv2")
GDSC1000 <- downloadPSet("GDSC1000")
CCLE <- downloadPSet("CCLE")
gCSI <- downloadPSet("gCSI")
PsetVec <- list(CCLE, gCSI, CTRPv2, GDSC1000)
names(PsetVec) <- c("CCLE", "gCSI",
                    "CTRPv2", "GDSC1000") 
############################## Combination of significance of interactions between datasets
Directory <- "~/Desktop/Drug_Tissue_Association/GSEA/"

if(Adjustment){
  CCLE_output <- readRDS(paste(Directory,"1_Adjusted_ResultList.rds", sep = ""))
  gCSI_output <- readRDS(paste(Directory,"2_Adjusted_ResultList.rds", sep = ""))
  CTRPv2_output <- readRDS(paste(Directory,"3_Adjusted_ResultList.rds", sep = ""))
  GDSC1000_output <- readRDS(paste(Directory,"4_Adjusted_ResultList.rds", sep = ""))
  Output_List <- list(CCLE_output, gCSI_output,
                      CTRPv2_output, GDSC1000)
}else{
  CCLE_output <- readRDS(paste(Directory, "1_ResultList.rds", sep = ""))
  gCSI_output <- readRDS(paste(Directory, "2_ResultList.rds", sep = ""))
  CTRPv2_output <- readRDS(paste(Directory, "3_ResultList.rds", sep = ""))
  GDSC1000 <- readRDS(paste(Directory, "4_ResultList.rds", sep = ""))
  Output_List <- list(CCLE_output, gCSI_output,
                      CTRPv2_output, GDSC1000)
}

names(Output_List) <- c("CCLE", "gCSI",
                        "CTRPv2", "GDSC1000")
###################
DrugTissueList <- list()
DrugTissue_Names <- c()
for(PsetName in names(PsetVec)){
  print(PsetName)
  
  TargetPSet <- PsetVec[[PsetName]]
  AUCmat <- summarizeSensitivityProfiles(TargetPSet,
                                         sensitivity.measure="auc_recomputed")
  DrugVec <- tolower(rownames(AUCmat))
  
  Pset_Result <- Output_List[[PsetName]]
  EnrichmentMat <- Pset_Result$Enrichment
  PvalMat <- Pset_Result$Pvalue
  TissueMat <- Pset_Result$Tissues
  
  TissueNum <- nrow(TissueMat)
  DrugTissueMat <- c()
  for(DrugIter in 1:nrow(AUCmat)){
    
    DrugTissueMat <- rbind(DrugTissueMat,
                           cbind(tolower(TissueMat[,DrugIter]),
                                 rep(DrugVec[DrugIter], TissueNum),
                                 PvalMat[,DrugIter],
                                 EnrichmentMat[,DrugIter]))
    DrugTissue_Names <- rbind(DrugTissue_Names, cbind(TissueMat[,DrugIter],
                                                      rep(DrugVec[DrugIter], TissueNum)))
  }
  colnames(DrugTissueMat) <- c("Tissue", "Drug",
                               "Pvalue", "Enrichment")
  DrugTissueList[[PsetName]] <- DrugTissueMat
}

colnames(DrugTissue_Names) <- c("Tissue", "Drug")
DrugTissue_Names <- DrugTissue_Names[-which(
  duplicated(paste(DrugTissue_Names[,"Tissue"],
                   DrugTissue_Names[,"Drug"], sep = " "))),]
######################
DrugTissue_PvalEnrich <- cbind(DrugTissue_Names,
                               matrix(rep(NA, nrow(DrugTissue_Names)*8),
                                      ncol = 8))
colnames(DrugTissue_PvalEnrich) <- c("Tissue", "Drug",
                                     paste("Pval", names(PsetVec), sep = "_"),
                                     paste("Enrichment", names(PsetVec),
                                           sep = "_"))
  
rownames(DrugTissue_PvalEnrich) <- paste(DrugTissue_PvalEnrich[,"Tissue"],
                                         DrugTissue_PvalEnrich[,"drug"],
                                         sep = "_")
########################
for(PsetName in names(PsetVec)){
  DrugTissueMat <- DrugTissueList[[PsetName]]
  
  MatchRows <- paste(DrugTissueMat[,"Tissue"],
                     DrugTissueMat[,"Drug"],
                     sep = "")
  
  DrugTissue_PvalEnrich[MatchRows,paste("Pval", PsetName, sep = "_")] <- DrugTissueMat[,"Pvalue"]
  DrugTissue_PvalEnrich[MatchRows,paste("Enrichment", PsetName, sep = "_")] <- DrugTissueMat[,"Enrichment"]
}

Combined_pval <- c()
for(DrugTissueIter in 1:nrow(DrugTissue_PvalEnrich)){
  Combined_pval <- c(Combined_pval,
                     combine.test(c(DrugTissue_PvalEnrich[,c(paste("Pval", names(PsetVec), sep = "_"))]),
                                  weight = rep(1,4), method = "z.transform"))
}
DrugTissue_PvalEnrich <- cbind(DrugTissue_PvalEnrich,
                               Combined_pval,
                               p.adjust(Combined_pval, method = "fdr"))
################
for(PsetName in names(PsetVec)){
  print(PsetName)
  DrugTissue_PvalEnrich[,paste("Pval", names(PsetVec), sep = "_")] <- p.adjust(DrugTissue_PvalEnrich[,paste("Pval", names(PsetVec), sep = "_")],
                                                                               method = "fdr")
}
################
colnames(DrugTissue_PvalEnrich) <- c("Tissue", "Drug",
                                     paste("FDR", names(PsetVec), sep = "_"),
                                     paste("Enrichment", names(PsetVec),
                                           "Combined_Pval", "Combined_FDR",sep = "_"))