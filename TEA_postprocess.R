########################
### functions
########################

postprocessingTEA <- function(ResultFileNames, PsetVec, Adjustment, GSEADir){
  
  ######################## Cell line number matrices
  CCLE_cclMat <- readRDS(file.path(GSEADir, "CCLE_originalAUC_Ncelline.rds"))
  gCSI_cclMat <- readRDS(file.path(GSEADir, "gCSI_originalAUC_Ncelline.rds"))
  CTRPv2_cclMat <- readRDS(file.path(GSEADir, "CTRPv2_originalAUC_Ncelline.rds"))
  GDSC1000_cclMat <- readRDS(file.path(GSEADir, "GDSC1000_originalAUC_Ncelline.rds"))
  
  CellNum_List <- list(CCLE_cclMat, gCSI_cclMat,
                       CTRPv2_cclMat, GDSC1000_cclMat)
  names(CellNum_List) <- names(PsetVec)
  ######################## Result lists
  CCLE_output <- readRDS(file.path(GSEADir,ResultFileNames["CCLE"]))
  gCSI_output <- readRDS(file.path(GSEADir,ResultFileNames["gCSI"]))
  CTRPv2_output <- readRDS(file.path(GSEADir,ResultFileNames["CTRPv2"]))
  GDSC1000_output <- readRDS(file.path(GSEADir,ResultFileNames["GDSC1000"]))
  Output_List <- list(CCLE_output, gCSI_output, CTRPv2_output, GDSC1000_output)
  
  names(Output_List) <- names(PsetVec)
  ################### Flattened matrix for each PSet
  DrugTissueList <- list()
  DrugTissue_Names <- c()
  for(PsetName in names(PsetVec)){
    print(PsetName)
    
    Pset_Result <- Output_List[[PsetName]]
    EnrichmentMat <- Pset_Result$Enrichment
    PvalMat <- Pset_Result$Pvalue
    
    CellLineMat <- CellNum_List[[PsetName]]
    
    DrugVec <- colnames(PvalMat)
    TissueVec <- rownames(PvalMat)
    
    TissueNum <- nrow(PvalMat)
    DrugTissueMat <- c()
    
    for(DrugIter in c(1:ncol(PvalMat))){
      
      NAind <- as.numeric(which(is.na(PvalMat[,DrugIter])))
      DrugTissueMat <- rbind(DrugTissueMat,
                             cbind(TissueVec[-NAind],
                                   rep(DrugVec[DrugIter], (TissueNum-length(NAind))),
                                   PvalMat[-NAind,DrugIter],
                                   EnrichmentMat[-NAind,DrugIter],
                                   CellLineMat[-NAind,DrugIter]))
      DrugTissue_Names <- rbind(DrugTissue_Names, cbind(TissueVec[-NAind],
                                                        rep(DrugVec[DrugIter], (TissueNum-length(NAind)))))
    }
    
    
    colnames(DrugTissueMat) <- c("Tissue", "Drug",
                                 "Pvalue", "Enrichment", "CCL_Num")
    DrugTissueList[[PsetName]] <- DrugTissueMat
  }
  
  colnames(DrugTissue_Names) <- c("Tissue", "Drug")
  DrugTissue_Names <- DrugTissue_Names[-which(
    duplicated(paste(DrugTissue_Names[,"Tissue"],
                     DrugTissue_Names[,"Drug"], sep = " "))),]
  ###################### Making the empty flattened matrix of combined PSets
  DrugTissue_PvalEnrich <- cbind(DrugTissue_Names,
                                 matrix(rep(NA, nrow(DrugTissue_Names)*12),
                                        ncol = 12))
  colnames(DrugTissue_PvalEnrich) <- c("Tissue", "Drug",
                                       paste("Pval", names(PsetVec), sep = "_"),
                                       paste("Enrichment", names(PsetVec),sep = "_"),
                                       paste("CCL_Num", names(PsetVec),sep = "_"))
  
  rownames(DrugTissue_PvalEnrich) <- paste(DrugTissue_PvalEnrich[,"Tissue"],
                                           DrugTissue_PvalEnrich[,"Drug"],
                                           sep = "_")
  ######################## Combining PSets in the flattened matrix
  for(PsetName in names(PsetVec)){
    DrugTissueMat <- DrugTissueList[[PsetName]]
    
    MatchRows <- paste(DrugTissueMat[,"Tissue"],
                       DrugTissueMat[,"Drug"],
                       sep = "_")
    
    DrugTissue_PvalEnrich[MatchRows,paste("Pval", PsetName, sep = "_")] <- DrugTissueMat[,"Pvalue"]
    DrugTissue_PvalEnrich[MatchRows,paste("Enrichment", PsetName, sep = "_")] <- DrugTissueMat[,"Enrichment"]
    DrugTissue_PvalEnrich[MatchRows,paste("CCL_Num", PsetName, sep = "_")] <- DrugTissueMat[,"CCL_Num"]
  }
  #################### combining Pvalues and FDR correction
  Combined_pval <- c()
  
  for(DrugTissueIter in 1:nrow(DrugTissue_PvalEnrich)){
    pvals <- as.numeric(c(DrugTissue_PvalEnrich[DrugTissueIter,
                                                c(paste("Pval", names(PsetVec),
                                                        sep = "_"))]))
    CCLNums <- as.numeric(c(DrugTissue_PvalEnrich[DrugTissueIter,
                                                  c(paste("CCL_Num", names(PsetVec),
                                                          sep = "_"))]))
    NotNa <- which(!is.na(pvals))
    if(length(NotNa) > 1){
      Combined_pval <- c(Combined_pval,
                         combine.test(pvals[NotNa],
                                      weight = CCLNums[NotNa],
                                      method = "z.transform"))
    }else{
      Combined_pval <- c(Combined_pval, NA)
    }
    
  }
  DrugTissue_PvalEnrich <- cbind(DrugTissue_PvalEnrich,
                                 Combined_pval, p.adjust(Combined_pval, method = "fdr"))
  ################ FDR correction for each PSet
  FDRMat <- c()
  for(PsetName in names(PsetVec)){
    print(PsetName)
    FDRMat <- cbind(FDRMat, p.adjust(DrugTissue_PvalEnrich[,paste("Pval", PsetName, sep = "_")],
                                     method = "fdr"))
  }
  
  DrugTissue_PvalEnrich <- cbind(DrugTissue_PvalEnrich[,c(1,2)],
                                 FDRMat, DrugTissue_PvalEnrich[,c(3:ncol(DrugTissue_PvalEnrich))])
  
  ################ Saving the files
  colnames(DrugTissue_PvalEnrich) <- c("Tissue", "Drug",
                                       paste("Pval", names(PsetVec), sep = "_"),
                                       paste("FDR", names(PsetVec), sep = "_"),
                                       paste("Enrichment", names(PsetVec), sep = "_"),
                                       paste("CCL_Num", names(PsetVec), sep = "_"),
                                       "Combined_Pval", "Combined_FDR")
  DrugTissue_PvalEnrich <- as.data.frame(DrugTissue_PvalEnrich)
  
  saveRDS(DrugTissue_PvalEnrich, file=file.path(GSEADir, sprintf("DrugTissue_PvalEnrich_%s.rds", ifelse(Adjustment, "adjustedAUC", "originalAUC"))))
  
  return(DrugTissue_PvalEnrich)
}

########################

### Original
ResultFileNames <- c("CCLE_originalAUC_ResultList.rds",
                     "gCSI_originalAUC_ResultList.rds",
                     "CTRPv2_originalAUC_ResultList.rds",
                     "GDSC1000_originalAUC_ResultList.rds")
names(ResultFileNames) <- c("CCLE", "gCSI", "CTRPv2", "GDSC1000")
Original_Results <- postprocessingTEA(ResultFileNames, PsetVec, Adjustment=FALSE, GSEADir)

### Adjusted
ResultFileNames <- c("CCLE_adjustedAUC_ResultList.rds",
                     "gCSI_adjustedAUC_ResultList.rds",
                     "CTRPv2_adjustedAUC_ResultList.rds",
                     "GDSC1000_adjustedAUC_ResultList.rds")
names(ResultFileNames) <- c("CCLE", "gCSI", "CTRPv2", "GDSC1000")
Adjusted_Results <- postprocessingTEA(ResultFileNames, PsetVec, Adjustment=TRUE, GSEADir)

tt <- list("Drug_Tissue_OriginalAUC"=Original_Results, "Drug_Tissue_AdjustedAUC"=Adjusted_Results)
WriteXLS::WriteXLS("tt", ExcelFileName=file.path(GSEADir, "DrugTissueAssocs_All.xlsx"))


### proportion of drugs and tissues involved in a significant association
tt <- Original_Results[!is.na(Original_Results[ , "Combined_FDR"]), ]
tt2 <- Adjusted_Results[!is.na(Original_Results[ , "Combined_FDR"]), ]
tt2 <- tt2[rownames(tt), ]
tt3 <- as.numeric(apply(cbind(tt[ , "Combined_FDR"], tt2[ , "Combined_FDR"]), 1, min, na.rm=TRUE))
names(tt3) <- rownames(tt)
iix <- which(tt3 < FDRcutoff)
length(sort(unique(tt[names(tt3)[iix], "Drug"]))) / length(sort(unique(tt[ , "Drug"])))

ncellines <- apply(data.matrix(tt[names(tt3), c("CCL_Num_CCLE", "CCL_Num_gCSI", "CCL_Num_CTRPv2", "CCL_Num_GDSC1000")]), 1, sum, na.rm=TRUE)
cor.test(ncellines, tt3, method="spearman", exact=FALSE)

### end

