
### list of drug-tissue associations and their statistics
drugTissueAssocs <- readRDS(file.path(GSEADir, "DrugTissue_PvalEnrich_Adjusted.rds"))
drugTissueAssocs <- data.frame(drugTissueAssocs)

### select the significant associations in meta-analysis
iix <- !is.na(drugTissueAssocs[ , "Combined_FDR"]) & drugTissueAssocs[ , "Combined_FDR"] < 0.05
drugTissueAssocs <- drugTissueAssocs[iix, , drop=FALSE]
drugTissueAssocs <- drugTissueAssocs[order(drugTissueAssocs[ , "Combined_FDR"], decreasing=FALSE), , drop=FALSE]

udrug <- sort(unique(drugTissueAssocs[ , "Drug"]))

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


ci <- NULL
for (i in 1:nrow(drugTissueAssocs)) {
  ### drug-tissue associations (d,t)
  dd <- drugTissueAssocs[i, "Drug"]
  tt <- drugTissueAssocs[i, "Tissue"]
  
  cis <- lapply(AUCList, function (x, drug, tissue, celltissue) {
    if (drug %in% rownames(x)) {
      ### create binary variable for (d,t)
      bb <- as.numeric(!is.na(celltissue[colnames(x)]) & celltissue[colnames(x)] == tissue)
      names(bb) <- colnames(x)
      aa <- x[drug, ]
      # (mRMRe::correlate(X=aa, Y=factor(bb, ordered=TRUE), method="cindex", outX=TRUE)$estimate / 2) + 0.5
      cindex <- survcomp::concordance.index(x=-bb, surv.time=aa, surv.event=rep(1, length(aa)), method="noether", na.rm=TRUE)
      res <- cindex[c(1, 2)]
    } else {
      res <- NULL
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











### end

