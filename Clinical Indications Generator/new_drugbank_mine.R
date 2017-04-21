source("get_drug_indications.R")

drug_indications<- get_drug_indications(unique(c(p_adj_false$Drug, p_adj_true$Drug)))

output2 <- melt(drug_indications, id.vars = c("drug"))
output2 <- na.omit(output2)
output2 <- output2[output2$variable != "Indication Browse", ]
output2 <- output2[output2$variable != "Liver Transplant Rejection",]
#output3 <- output2[output2$variable %in% unique(c(p_adj_false$Tissue, p_adj_true$Tissue)) | grepl("oma", output2$variable, ignore.case = TRUE) | grepl("cancer", output2$variable, ignore.case = TRUE), ]

output3$variable <- as.character(output3$variable)


output3[grep("renal", output3$variable, ignore.case = TRUE), "variable"] <- "kidney"
output3[grep("breast", output3$variable, ignore.case = TRUE), "variable"] <- "breast"
output3[grep("thyroid", output3$variable, ignore.case = TRUE), "variable"] <- "thyroid"
output3[grep("pancreatic", output3$variable, ignore.case = TRUE), "variable"] <- "pancreas"
output3[grep("ovarian", output3$variable, ignore.case = TRUE), "variable"] <- "ovary"
output3[grep("prostate", output3$variable, ignore.case = TRUE), "variable"] <- "prostate"

output3[grep("lymphoma", output3$variable, ignore.case = TRUE), "variable"] <- "haematopoietic_and_lymphoid_tissue"
output3[grep("myeloma", output3$variable, ignore.case = TRUE), "variable"] <- "haematopoietic_and_lymphoid_tissue"
output3[grep("thymoma", output3$variable, ignore.case = TRUE), "variable"] <-  "haematopoietic_and_lymphoid_tissue"
output3[grep("thymic", output3$variable, ignore.case = TRUE), "variable"] <-  "haematopoietic_and_lymphoid_tissue"
output3[grep("leukemia", output3$variable, ignore.case = TRUE), "variable"] <- "haematopoietic_and_lymphoid_tissue"
output3[grep("leukaemia", output3$variable, ignore.case = TRUE), "variable"] <- "haematopoietic_and_lymphoid_tissue"

output3[grep("bone", output3$variable, ignore.case = TRUE), "variable"] <- "bone"
output3[grep("osteo", output3$variable, ignore.case = TRUE), "variable"] <- "bone"

output3[grep("head and neck", output3$variable, ignore.case = TRUE), "variable"] <- "head_and_neck"
output3[grep("cervical", output3$variable, ignore.case = TRUE), "variable"] <- "cervix"

output3[grep("sarcoma", output3$variable, ignore.case = TRUE), "variable"] <-  "soft_tissue"
output3[grep("rhabdomyosarcoma", output3$variable, ignore.case = TRUE), "variable"] <- "soft_tissue"
output3[grep("kaposi", output3$variable, ignore.case = TRUE), "variable"] <- "soft_tissue"
output3[grep("soft tissue", output3$variable, ignore.case = TRUE), "variable"] <- "soft_tissue"

output3[grep("bladder", output3$variable, ignore.case = TRUE), "variable"] <- "urinary_tract"

output3[grep("uterine", output3$variable, ignore.case = TRUE), "variable"] <- "endometrium"
output3[grep("choriocarcinoma", output3$variable, ignore.case = TRUE), "variable"] <- "endometrium"

output3[grep("non small cell lung", output3$variable, ignore.case = TRUE), "variable"] <- "NSCLC"
output3[grep("non-small cell lung", output3$variable, ignore.case = TRUE), "variable"] <- "NSCLC"
output3[grep("non-small-cell lung", output3$variable, ignore.case = TRUE), "variable"] <- "NSCLC"
output3[grep("squamous cell lung", output3$variable, ignore.case = TRUE), "variable"] <- "NSCLC"
output3[grep("small cell lung", output3$variable, ignore.case = TRUE), "variable"] <- "SCLC"

output3[grep("neuro", output3$variable, ignore.case = TRUE), "variable"] <- "central_nervous_system"
output3[grep("chordomas", output3$variable, ignore.case = TRUE), "variable"] <- "central_nervous_system"
output3[grep("meningeal", output3$variable, ignore.case = TRUE), "variable"] <-  "central_nervous_system"
output3[grep("glioblastoma", output3$variable, ignore.case = TRUE), "variable"] <-  "central_nervous_system"

output3[grep("gastrointestinal", output3$variable, ignore.case = TRUE), "variable"] <- "large_intestine"
output3[grep("gastric", output3$variable, ignore.case = TRUE), "variable"] <- "large_intestine"
output3[grep("colorectal", output3$variable, ignore.case = TRUE), "variable"] <-  "large_intestine"

output3[grep("merkel", output3$variable, ignore.case = TRUE), "variable"] <-  "skin"
output3[grep("melanoma", output3$variable, ignore.case = TRUE), "variable"] <- "skin"

output3[grep("hepatocellular", output3$variable, ignore.case = TRUE), "variable"] <-  "liver"

output3 <- output3[output3$variable %in% unique(c(p_adj_false$Tissue, p_adj_true$Tissue)), ]

write.xlsx(unique(output3)[, c("drug", "variable")], file = "new_clinical_indications.xlsx")

saveRDS(output3, file = "new_clinical_indications.RDS")

all_drugs <- get_drug_indications(unique(c(colnames(one_offs), colnames(combined1))))

"727 drugs"
"91 drugs in drugbank"

all_drugs_melt <- melt(all_drugs, id.vars = c("drug"))
cancer_drugs <- c()
not_cancer_drugs <- c()

for(x in unique(all_drugs_melt$drug)){
  set <- all_drugs_melt[all_drugs_melt$drug == x, ]
  set <- na.omit(set)
  set$variable <- tolower(set$variable)
  if(TRUE %in% grepl("oma", set$variable) || TRUE %in% grepl("cancer", set$variable))
  {
    cancer_drugs <- c(cancer_drugs, x)
  } else { not_cancer_drugs <- c(not_cancer_drugs, x)}
}

#44 cancer drugs
#47 not cancer drugs 

