pres <- data.frame(matrix(ncol = 5, nrow = ncol(combined1)))
rownames(pres) <- colnames(combined1)
colnames(pres) <- c("CTRPv2", "CCLE", "GDSC", "gCSI", "NCI60")

pres[, "CTRPv2"] <- rownames(pres) %in% rownames(CTRPv2@drug)
pres[, "CCLE"] <- rownames(pres) %in% rownames(CCLE@drug)
pres[, "GDSC"] <- rownames(pres) %in% rownames(GDSC@drug)
pres[, "gCSI"] <- rownames(pres) %in% rownames(gCSI@drug)
pres[, "NCI60"] <- rownames(pres) %in% NCI60@curation$drug$unique.drugid

write.csv(pres, file = "suppfile2.csv")