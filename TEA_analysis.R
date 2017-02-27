### tissue enrichment analysis

message("-------------------------\nTissue enrichment analysis\n-------------------------")

for(PsetIter in 1:length(PsetVec)){

  message(sprintf("%s PharmacoSet [%i drugs]", names(PsetVec[PsetIter]), length(drugNames(PsetVec[[PsetIter]]))))
  ### file to store each results
  resfn <- file.path(GSEADir, sprintf("%s_%s", names(PsetVec)[PsetIter], ifelse(Adjustment, "adjustedAUC", "originalAUC")))
  if (!file.exists(paste(resfn, "ResultList.rds", sep="_"))) {
  
    TargetPSet <- PsetVec[[PsetIter]]
    AUCmat <- PharmacoGx::summarizeSensitivityProfiles(TargetPSet, sensitivity.measure="auc_recomputed")
    
    if(Adjustment) {
      AUCFilled <- Hmisc::impute(AUCmat, fun=median)
      class(AUCFilled) <- class(AUCmat)
    }

    DrugVec <- rownames(AUCmat)
    tissue <- as.character(cellInfo(TargetPSet)[ ,"tissueid_TEA"])
    utissue <- sort(unique(tissue))
  
    pb <- utils::txtProgressBar(min=0, max=length(DrugVec), style=3)
  
    AUCmat.adj <- matrix(NA, nrow=nrow(AUCmat), ncol=ncol(AUCmat), dimnames=dimnames(AUCmat))
    
    EnrichmentMat <- PvalueMat <- NcellineMat <- matrix(NA, ncol=length(DrugVec), nrow=length(utissue), dimnames=list(utissue, DrugVec))
    
    for(DrugIter in 1:length(DrugVec)) {
    
      AUCmat_tissueinf <- data.frame(AUCmat[DrugVec[DrugIter], ], tissue)
      colnames(AUCmat_tissueinf) <- c("AUC", "tissueid")
      if (sum(complete.cases(AUCmat_tissueinf)) >= 10) {
        ### PCA adjustment
        if(Adjustment) {
          CorVec <- c()
          for(DrugIter2 in 1:nrow(AUCmat)) {
            if (sum(complete.cases(AUCmat[DrugVec[DrugIter], ], AUCmat[DrugIter2, ])) >= 10) {
              CorVec <- c(CorVec, cor(AUCmat[DrugVec[DrugIter], ], AUCmat[DrugIter2, ], method = "spearman", use="complete.obs"))
            } else {
              CorVec <- c(CorVec, NA)
            }
          }
          LowCorDrugs <- which(CorVec < max(quantile(na.omit(CorVec))[2], 0))
          PCA <- prcomp(t(AUCFilled[LowCorDrugs, ]))
          AUCmat_tissueinf[ , "AUC"] <- AUCmat_tissueinf[ , "AUC"] - as.numeric(PCA$x[ , 1])
          AUCmat.adj[DrugVec[DrugIter], ] <- AUCmat[DrugVec[DrugIter], ] - as.numeric(PCA$x[ , 1])
        }
      
        AUCmat_tissueinf[,"AUC"] <- (AUCmat_tissueinf[,"AUC"] - median(AUCmat_tissueinf[,"AUC"], na.rm=TRUE)) / mad(AUCmat_tissueinf[,"AUC"], na.rm=TRUE)
        ### remove hematopoetic and lymphoid tissue
        AUCmat_tissueinf[which(grepl("lymphoid", AUCmat_tissueinf[ , "tissueid"])), "tissueid"] <- NA
        AUCmat_tissueinf <- na.omit(AUCmat_tissueinf)
        celline_tissueinf <- cbind("g"=rownames(AUCmat_tissueinf), "s"=AUCmat_tissueinf[ , "tissueid"])

        ### number of cell lines per tissue type
        tt <- table(celline_tissueinf[ , "s"])
        NcellineMat[names(tt), DrugVec[DrugIter]] <- tt

        genelevelstats <- AUCmat_tissueinf[ , "AUC", drop=FALSE]
        gsc1 <- piano::loadGSC(celline_tissueinf)

        if(sum(table(AUCmat_tissueinf[ , "tissueid"]) >= TissueSize[1] & table(AUCmat_tissueinf[ , "tissueid"]) <= TissueSize[2], na.rm=TRUE) > 1) {

          gsea_out <- piano::runGSA(geneLevelStats=genelevelstats, geneSetStat="gsea", gsc=gsc1,  nPerm=nperm + nbcore - (nperm %% nbcore), ncpus=nbcore, gsSizeLim=TissueSize, adjMethod="none", verbose=FALSE)

          gseares <- try(piano::GSAsummaryTable(gsea_out))
          if (class(gseares) != "try-error" && nrow(gseares) > 1) {

            ### get p-values and enrichment scores
            PvalueVec <- EnrichmentVec <- rep(NA, nrow(gseares))
            names(PvalueVec) <- names(EnrichmentVec) <- as.character(gseares[ , "Name"])
            ### p-values
            if ("p (dist.dir.up)" %in% colnames(gseares)) {
              iix <- is.na(PvalueVec) & !is.na(gseares[ , "p (dist.dir.up)"])
              PvalueVec[iix] <- gseares[iix, "p (dist.dir.up)"]
            }
            if ("p (dist.dir.dn)" %in% colnames(gseares)) {
              iix <- is.na(PvalueVec) & !is.na(gseares[ , "p (dist.dir.dn)"])
              PvalueVec[iix] <- 1
            }
            PvalueVec[which(PvalueVec == 0)] <- 1 / (nperm + nbcore - (nperm %% nbcore) + 1)
            if (length(PvalueVec) < 1 | is.null(PvalueVec)) {
              stop(sprintf("Error for drug %s", DrugVec[DrugIter]))
            }
            PvalueMat[names(PvalueVec), DrugVec[DrugIter]] <- PvalueVec
            ### enrichment scores
            if ("Stat (dist.dir)" %in% colnames(gseares)) {
              iix <- !is.na(gseares[ , "Stat (dist.dir)"])
              EnrichmentVec[iix] <- gseares[iix, "Stat (dist.dir)"]
            }
            EnrichmentMat[names(EnrichmentVec), DrugVec[DrugIter]] <- EnrichmentVec

          }
        }
      } 
      utils::setTxtProgressBar(pb, DrugIter)
    }
    message("")
    
    if (Adjustment) { AUCm <- AUCmat.adj } else { AUCm <- AUCmat }
    ResultList <- list("Enrichment"=EnrichmentMat, "Pvalue"=PvalueMat, "N"=NcellineMat, "AUC"=AUCm)
    saveRDS(object=ResultList, file=paste(resfn, "ResultList.rds", sep="_"))
    saveRDS(object=AUCm, file=paste(resfn, "AUC.rds", sep="_"))
    saveRDS(object=NcellineMat, file=paste(resfn, "Ncelline.rds", sep="_"))PsetIter <- 4
    rm(list=c("ResultList", "AUCmat.adj", "AUCmat", "AUCm", "EnrichmentMat", "PvalueMat", "NcellineMat"))
    gc()
  }
}

### end



