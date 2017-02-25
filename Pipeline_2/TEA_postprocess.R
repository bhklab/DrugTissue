########################
### functions
########################
postprocessingTEA <- function(Ncellline_Files, PsetVec, Adjustment, GSEADir) {
  CCLE_cclMat <- readRDS(paste(GSEADir, as.character(Ncellline_Files["CCLE"]), sep = "")) ## "CCLE_originalAUC_Ncelline.rds"
  gCSI_cclMat <- readRDS(paste(GSEADir, as.character(Ncellline_Files["gCSI"]), sep = "")) ## "gCSI_originalAUC_Ncelline.rds"
  CTRPv2_cclMat <- readRDS(paste(GSEADir, as.character(Ncellline_Files["CTRPv2"]), sep = "")) ## "CTRPv2_originalAUC_Ncelline.rds"
  GDSC1000_cclMat <- readRDS(paste(GSEADir, as.character(Ncellline_Files["GDSC1000"]), sep = "")) ## "GDSC1000_originalAUC_Ncelline.rds"
  
  CellNum_List <- list(CCLE_cclMat, gCSI_cclMat,
                       CTRPv2_cclMat, GDSC1000_cclMat)
  names(CellNum_List) <- names(PsetVec)
  ########################
  
  if(Adjustment){
    CCLE_output <- readRDS(paste(GSEADir,"CCLE_adjustedAUC_ResultList.rds", sep = ""))
    gCSI_output <- readRDS(paste(GSEADir,"gCSI_adjustedAUC_ResultList.rds", sep = ""))
    CTRPv2_output <- readRDS(paste(GSEADir,"CTRPv2_adjustedAUC_ResultList.rds", sep = ""))
    GDSC1000_output <- readRDS(paste(GSEADir,"GDSC1000_adjustedAUC_ResultList.rds", sep = ""))
    Output_List <- list(CCLE_output, gCSI_output,
                        CTRPv2_output, GDSC1000_output)
  }else{
    CCLE_output <- readRDS(paste(GSEADir, "CCLE_originalAUC_ResultList.rds", sep = ""))
    gCSI_output <- readRDS(paste(GSEADir, "gCSI_originalAUC_ResultList.rds", sep = ""))
    CTRPv2_output <- readRDS(paste(GSEADir, "CTRPv2_originalAUC_ResultList.rds", sep = ""))
    GDSC1000_output <- readRDS(paste(GSEADir, "GDSC1000_originalAUC_ResultList.rds", sep = ""))
    Output_List <- list(CCLE_output, gCSI_output,
                        CTRPv2_output, GDSC1000_output)
  }
  
  names(Output_List) <- names(PsetVec)
  ###################
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
  ######################
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
  ########################
  for(PsetName in names(PsetVec)){
    DrugTissueMat <- DrugTissueList[[PsetName]]
    
    MatchRows <- paste(DrugTissueMat[,"Tissue"],
                       DrugTissueMat[,"Drug"],
                       sep = "_")
    
    DrugTissue_PvalEnrich[MatchRows,paste("Pval", PsetName, sep = "_")] <- DrugTissueMat[,"Pvalue"]
    DrugTissue_PvalEnrich[MatchRows,paste("Enrichment", PsetName, sep = "_")] <- DrugTissueMat[,"Enrichment"]
    DrugTissue_PvalEnrich[MatchRows,paste("CCL_Num", PsetName, sep = "_")] <- DrugTissueMat[,"CCL_Num"]
  }
  ####################
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
  ################
  FDRMat <- c()
  for(PsetName in names(PsetVec)){
    print(PsetName)
    FDRMat <- cbind(FDRMat, p.adjust(DrugTissue_PvalEnrich[,paste("Pval", PsetName, sep = "_")],
                                     method = "fdr"))
  }
  
  DrugTissue_PvalEnrich <- cbind(DrugTissue_PvalEnrich[,c(1,2)],
                                 FDRMat, DrugTissue_PvalEnrich[,c(3:ncol(DrugTissue_PvalEnrich))])
  
  ################
  colnames(DrugTissue_PvalEnrich) <- c("Tissue", "Drug",
                                       paste("Pval", names(PsetVec), sep = "_"),
                                       paste("FDR", names(PsetVec), sep = "_"),
                                       paste("Enrichment", names(PsetVec), sep = "_"),
                                       paste("CCL_Num", names(PsetVec), sep = "_"),
                                       "Combined_Pval", "Combined_FDR")
  DrugTissue_PvalEnrich <- as.data.frame(DrugTissue_PvalEnrich)
  
  if(Adjustment){
    FileName <- "adjusted"
  }else{
    FileName <- "original"
  }
  saveRDS(DrugTissue_PvalEnrich, file = paste(GSEADir, "DrugTissue_PvalEnrich_",FileName,"AUC.rds"))
  WriteXLS(DrugTissue_PvalEnrich, file = paste(GSEADir, "DrugTissue_PvalEnrich_",FileName,"AUC.xls"))
  
}