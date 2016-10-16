pres <- data.frame(matrix(ncol = 5, nrow = ncol(combined1)))
rownames(pres) <- colnames(combined1)
colnames(pres) <- c("CTRPv2", "CCLE", "GDSC1000", "gCSI", "NCI60")

pres[, "CTRPv2"] <- rownames(pres) %in% rownames(CTRPv2@drug)
pres[, "CCLE"] <- rownames(pres) %in% rownames(CCLE@drug)
pres[, "GDSC1000"] <- rownames(pres) %in% rownames(GDSC1000@drug)
pres[, "gCSI"] <- rownames(pres) %in% rownames(gCSI@drug)
pres[, "NCI60"] <- rownames(pres) %in% NCI60@curation$drug$unique.drugid

pres[pres == TRUE] <- 1
pres[pres == FALSE] <- ""

write.xlsx(pres, file = "Supplementary_file_1.xlsx", rowNames = TRUE)
