
############################## Combination of significance of interactions between datasets
if(Adjustment){
  CCLE_output <- readRDS(paste(GSEADir,"1_Adjusted_ResultList.rds", sep = ""))
  gCSI_output <- readRDS(paste(GSEADir,"2_Adjusted_ResultList.rds", sep = ""))
  CTRPv2_output <- readRDS(paste(GSEADir,"3_Adjusted_ResultList.rds", sep = ""))
  GDSC1000_output <- readRDS(paste(GSEADir,"4_Adjusted_ResultList.rds", sep = ""))
  Output_List <- list(CCLE_output, gCSI_output,
                      CTRPv2_output, GDSC1000_output)
}else{
  CCLE_output <- readRDS(paste(GSEADir, "1_ResultList.rds", sep = ""))
  gCSI_output <- readRDS(paste(GSEADir, "2_ResultList.rds", sep = ""))
  CTRPv2_output <- readRDS(paste(GSEADir, "3_ResultList.rds", sep = ""))
  GDSC1000 <- readRDS(paste(GSEADir, "4_ResultList.rds", sep = ""))
  Output_List <- list(CCLE_output, gCSI_output,
                      CTRPv2_output, GDSC1000_output)
}

names(Output_List) <- c("CCLE", "gCSI",
                        "CTRPv2", "GDSC1000")
###################
DrugTissueList <- list()
DrugTissue_Names <- c()
for(PsetName in names(PsetVec)[3]){
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
                     survcomp::combine.test(c(DrugTissue_PvalEnrich[,c(paste("Pval", names(PsetVec), sep = "_"))]),
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