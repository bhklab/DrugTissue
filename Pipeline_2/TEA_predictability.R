########################
### functions
########################

predictabilityTEA <- function (drugTissueAssocs, FDRcutoff, PsetVec, Adjustment, quantileAUC, GSEADir) {
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
          qAUC <- quantile(aa[which(bb == 1)], probs=quantileAUC, na.rm=TRUE)
          cindex <- Hmisc::rcorr.cens(x=bb, S=Surv(aa, rep(1, length(aa))), outx=TRUE)
          res <- list("cindex"=cindex["C Index"], "se"=cindex["S.D."], "QuantileAUC"=qAUC)
        }
      }
      return (res)
    }, drug=dd, tissue=tt, celltissue=CellTissue)
    cis <- cis[!sapply(cis, is.null)]
    if (length(cis) < 1) {
      ci2 <- NA
    } else {
      ci2 <- survcomp::combine.est(x=sapply(cis, function (x) { return (x[[1]]) }), x.se=sapply(cis, function (x) { return (x[[2]]) }), hetero=TRUE)$estimate
    }
    ci <- rbind(ci, c("Cindex"=ci2, "quantileAUC"=median(sapply(cis, function (x) { return (x[[3]]) }), na.rm=TRUE)))
  }
  rownames(ci) <- rownames(drugTissueAssocs)
 
  return (ci)
}

########################

### list of drug-tissue associations and their statistics
### originalAUC
drugTissueAssocs.original <- readRDS(file.path(GSEADir, sprintf("DrugTissue_PvalEnrich_%s.rds", "originalAUC")))
drugTissueAssocs.original <- data.frame(drugTissueAssocs.original)
ci.original <- predictabilityTEA(drugTissueAssocs.original, FDRcutoff, PsetVec, Adjustment=FALSE, quantileAUC, GSEADir)
### adjustedAUC
drugTissueAssocs.adjusted <- readRDS(file.path(GSEADir, sprintf("DrugTissue_PvalEnrich_%s.rds", "adjustedAUC")))
drugTissueAssocs.adjusted <- data.frame(drugTissueAssocs.adjusted)
ci.adjusted <- predictabilityTEA(drugTissueAssocs.adjusted, FDRcutoff, PsetVec, Adjustment=TRUE, quantileAUC, GSEADir)


intersect(rownames(ci.original), rownames(ci.adjusted))

drugTissueAssocs <- cbind(drugTissueAssocs.original[rownames(ci.original), c("Tissue", "Drug", "Combined_Pval", "Combined_FDR")], ci.original, "Num_Cell_Lines"=apply(data.matrix(drugTissueAssocs.original[rownames(ci.original), c("CCL_Num_CCLE", "CCL_Num_gCSI", "CCL_Num_CTRPv2", "CCL_Num_GDSC1000")]), 1, sum, na.rm=TRUE), "Adjustment"="NO")
drugTissueAssocs <- rbind(drugTissueAssocs, cbind(drugTissueAssocs.adjusted[rownames(ci.adjusted), c("Tissue", "Drug", "Combined_Pval", "Combined_FDR")], ci.adjusted, "Num_Cell_Lines"=apply(data.matrix(drugTissueAssocs.adjusted[rownames(ci.adjusted), c("CCL_Num_CCLE", "CCL_Num_gCSI", "CCL_Num_CTRPv2", "CCL_Num_GDSC1000")]), 1, sum, na.rm=TRUE), "Adjustment"="YES"))

drugTissueAssocs <- drugTissueAssocs[order(rownames(drugTissueAssocs)), , drop=FALSE]


ll <- list("Drug Tissue Associations"=drugTissueAssocs)
WriteXLS::WriteXLS("ll", ExcelFileName=file.path(GSEADir, sprintf("drugTissueAssocs_FDR_%i.xlsx", ceiling(FDRcutoff*100))))










### end

