### tissue enrichment analysis

message("Tissue enrichment analysis")

for(PsetIter in 1:length(PsetVec)){

  message(sprintf("%s PharmacoSet [%i drugs]", names(PsetVec[PsetIter]), length(drugNames(PsetVec[[PsetIter]]))))
  ### file to store each results
  resfn <- file.path(GSEADir, sprintf("%s_%s_ResultsList.rds", names(PsetVec)[PsetIter], ifelse(Adjustment, "adjustedAUC", "originalAUC")))
  if (!file.exists(resfn)) {
  
    TargetPSet <- PsetVec[[PsetIter]]
    AUCmat <- PharmacoGx::summarizeSensitivityProfiles(TargetPSet, sensitivity.measure="auc_recomputed")

    if(Adjustment) {
      AUCFilled <- Hmisc::impute(AUCmat, fun=median)
      class(AUCFilled) <- class(AUCmat)
    }

    DrugVec <- rownames(AUCmat)
  
    pb <- utils::txtProgressBar(min=0, max=length(DrugVec), style=3)
    
    EnrichmentMat <- c()
    PvalMat <- c()
    TissueMat <- c()
  
    for(DrugIter in 1:length(DrugVec)) {
    
      AUCmat_tissueinf <- data.frame(AUCmat[DrugVec[DrugIter], ], cellInfo(TargetPSet)[ ,"tissueid_TEA"])
      colnames(AUCmat_tissueinf) <- c("AUC", "tissueid")
    
    
      ### PCA adjustment
      if(Adjustment) {
        CorVec <- c()
        for(DrugIter2 in 1:nrow(AUCmat)) {
          if (sum(complete.cases(AUCmat[DrugVec[DrugIter], ], AUCmat[DrugIter2, ])) > 10) {
            CorVec <- c(CorVec, cor(AUCmat[DrugVec[DrugIter], ], AUCmat[DrugIter2, ], method = "spearman", use="complete.obs"))
          } else {
            CorVec <- c(CorVec, NA)
          }
        }
        LowCorDrugs <- which(CorVec < max(quantile(na.omit(CorVec))[2], 0))
        PCA <- prcomp(t(AUCFilled[LowCorDrugs, ]))
        AUCmat_tissueinf[ , "AUC"] <- (AUCmat_tissueinf[ , "AUC"] - as.numeric(PCA$x[ , 1]))
      }
      
      AUCmat_tissueinf[,"AUC"] <- (AUCmat_tissueinf[,"AUC"] - median(AUCmat_tissueinf[,"AUC"], na.rm=TRUE)) / mad(AUCmat_tissueinf[,"AUC"], na.rm=TRUE)
      ### remove hematopoetic and lymphoid tissue
      AUCmat_tissueinf[which(grepl("lymphoid", AUCmat_tissueinf[ , "tissueid"])), "tissueid"] <- NA
      AUCmat_tissueinf <- na.omit(AUCmat_tissueinf)
      celline_tissueinf <- cbind("g"=rownames(AUCmat_tissueinf), "s"=AUCmat_tissueinf[ , "tissueid"])
    
      genelevelstats <- AUCmat_tissueinf[ , "AUC", drop=FALSE]
      gsc1 <- piano:::loadGSC(celline_tissueinf)
    
      if(sum(table(AUCmat_tissueinf[ , "tissueid"]) >= TissueSize[1] & table(AUCmat_tissueinf[ , "tissueid"]) <= TissueSize[2], na.rm=TRUE) > 1) {
       
        # gsea_out <- piano::runGSA(geneLevelStats=genelevelstats, geneSetStat="gsea", gsc=gsc1,  nPerm=nperm + nbcore - (nperm %% nbcore), ncpus=nbcore, gsSizeLim=TissueSize, adjMethod="none", verbose=FALSE)
      
        gseares <- piano::GSAsummaryTable(gsea_out)
        gseares <- cbind(gseares, "p"=NA, "p adj"=NA)
        gseares$`p (dist.dir.up)`
        PvalueVec <- gseares$`p (dist.dir.up)`
        NAind <- which(is.na(PvalueVec) & !is.na(gseares$`p (dist.dir.dn)`))
      
        PvalueVec[NAind] <- 1
        PvalueVec[which(PvalueVec == 0)] <- 1/(nperm + nbcore - (nperm %% nbcore) + 1)
      
        EnrichmentVec <- gseares$`Stat (dist.dir)`
      
        TissueVec <- gseares$Name
      
        EnrichmentMat <- cbind(EnrichmentMat, EnrichmentVec)
        PvalMat <- cbind(PvalMat, PvalueVec)
        TissueMat <- cbind(TissueMat, TissueVec)
      } else {
        EnrichmentMat <- cbind(EnrichmentMat, NA)
        PvalMat <- cbind(PvalMat, NA)
        TissueMat <- cbind(TissueMat, NA)
      }
      
      utils::setTxtProgressBar(pb, DrugIter)
    }
    message("\n")

    ResultList <- list(EnrichmentMat, PvalMat, TissueMat)
    names(ResultList) <- c("Enrichment", "Pvalue", "Tissues")
    saveRDS(object=ResultList, file=myfn)
    rm(ResultList)
    gc()
    
  }
  
}

### end



