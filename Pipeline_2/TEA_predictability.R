########################
### functions
########################

predictabilityTEA <- function (drugTissueAssocs, FDRcutoff, PsetVec, Adjustment, GSEADir) {
  ### select the significant associations in meta-analysis
  iix <- !is.na(drugTissueAssocs[ , "Combined_FDR"]) & drugTissueAssocs[ , "Combined_FDR"] < FDRcutoff
  drugTissueAssocs <- drugTissueAssocs[iix, , drop=FALSE]
  drugTissueAssocs <- drugTissueAssocs[order(drugTissueAssocs[ , "Combined_FDR"], decreasing=FALSE), , drop=FALSE]

  udrug <- sort(unique(drugTissueAssocs[ , "Drug"]))

  ### estract AUC values
  AUCList <- CellTissue <- NULL
  for(PsetIter in 1:length(PsetVec)) {
    resfn <- file.path(GSEADir, sprintf("%s_%s", names(PsetVec)[PsetIter], ifelse(Adjustment, "adjustedAUC", "originalAUC")))
    AUCm <- readRDS(file=paste(resfn, "AUC.rds", sep="_"))
    AUCList <- c(AUCList, list(AUCm))
    tt <- PharmacoGx::cellInfo(PsetVec[[PsetIter]])[ , "tissueid_TEA", drop=FALSE]
    tt[which(grepl("lymphoid", tt[ , 1])), 1] <- NA
    tt2 <- rownames(tt)
    names(tt2) <- as.character(tt[ , 1])
    CellTissue <- c(CellTissue, tt2)
  }
  names(AUCList) <- names(PsetVec)
  CellTissue <- CellTissue[!duplicated(CellTissue)]
  tt <- names(CellTissue)
  names(tt) <- CellTissue
  CellTissue <- tt

  ### compute concordance indices
  ci <- NULL
  for (i in 1:nrow(drugTissueAssocs)) {
    ### drug-tissue associations (d,t)
    dd <- drugTissueAssocs[i, "Drug"]
    tt <- drugTissueAssocs[i, "Tissue"]
  
    cis <- lapply(AUCList, function (x, drug, tissue, celltissue) {
      res <- NULL
      if (drug %in% rownames(x)) {
        ### create binary variable for (d,t)
        bb <- as.numeric(!is.na(celltissue[colnames(x)]) & celltissue[colnames(x)] == tissue)
        names(bb) <- colnames(x)
        aa <- x[drug, ]
        if (sum(complete.cases(bb, aa)) >= 10) {
          cindex <- Hmisc::rcorr.cens(x=bb, S=Surv(aa, rep(1, length(aa))), outx=TRUE)
          res <- list("cindex"=cindex["C Index"], "se"=cindex["S.D."])
        }
      }
      return (res)
    }, drug=dd, tissue=tt, celltissue=CellTissue)
    cis <- cis[!sapply(cis, is.null)]
    if (length(cis) < 1) {
      cis <- NA
    } else {
      cis <- survcomp::combine.est(x=sapply(cis, function (x) { return (x[[1]]) }), x.se=sapply(cis, function (x) { return (x[[2]]) }), hetero=TRUE)$estimate
    }
    ci <- c(ci, cis)
  }
  names(ci) <- rownames(drugTissueAssocs)
 
  return (ci)
}

########################

### list of drug-tissue associations and their statistics
drugTissueAssocs <- readRDS(file.path(GSEADir, sprintf("DrugTissue_PvalEnrich_%s.rds", ifelse(Adjustment, "adjustedAUC", "originalAUC"))))
drugTissueAssocs <- data.frame(drugTissueAssocs)

predictabilityTEA(drugTissueAssocs, FDRcutoff, PsetVec, Adjustment, GSEADir)



drugTissueAssocs <- cbind(drugTissueAssocs, "tissue_cindex"=ci)
ll <- list("")
WriteXLS::WriteXLS










### end

