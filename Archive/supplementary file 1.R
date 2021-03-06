.libPaths("/rlibs/")
library(PharmacoGx)

if(!file.exists("./temp/combined1.RData"))
{
  source("./gsea_with_AUC.R")
} else {
  load("./temp/combined1.RData")
}

if(!file.exists("./temp/amla.RData"))
{
  source("./gsea_with_AUC.R")
} else {
  load("./temp/amla.RData")
}

CCLE <- downloadPSet("CCLE_2013")
GDSC100 <- downloadPSet("GDSC1000")
gCSI <- downloadPSet("gCSI")
CTRPv2 <- downloadPSet("CTRPv2")


pres <- data.frame(matrix(ncol = 5, nrow = ncol(combined1)))
rownames(pres) <- colnames(combined1)
colnames(pres) <- c("CTRPv2", "CCLE", "GDSC1000", "gCSI", "NCI60")

pres[, "CTRPv2"] <- rownames(pres) %in% colnames(ml$CTRPv2)
pres[, "CCLE"] <- rownames(pres) %in% colnames(ml$CCLE)
pres[, "GDSC1000"] <- rownames(pres) %in% colnames(ml$GDSC1000)
pres[, "gCSI"] <- rownames(pres) %in% colnames(ml$gCSI)
pres[, "NCI60"] <- rownames(pres) %in% NCI60@curation$drug$unique.drugid

pres[pres == TRUE] <- 1
pres[pres == FALSE] <- ""
pres["Afatinib", "NCI60"] <- 1

write.xlsx(pres, file = "./output/Supplementary_file_1.xlsx", row.names = TRUE, col.names = TRUE)
