wordmine <- data.frame()

inputdruglist <- unique(c(rownames(gCSI@drug), rownames(CCLE@drug), rownames(GDSC1000@drug), rownames(CTRPv2@drug), breast, colorectal, prostate))
inputdruglist <- unique(toupper(gsub(badchars, "", inputdruglist)))

for(x in inputdruglist)
{
  if(! x %in% wordmine[1,])
  {
    out <- getDrugTests(x)
  }
  
  if(!is.na(out))
  {
    for(y in 1:nrow(out))
    {
      wordmine[x, out[y, "X2"]] <- out[y, "X1"]
    }
  }
}

wordmine <- as.data.frame(wordmine)
wordmine[is.na(wordmine)] <- 0
wordmine <- wordmine[colnames(wordmine) != "other"]
wordmine <- wordmine[colnames(wordmine) != "large"]
wordmine <- wordmine[colnames(wordmine) != "tongue"]
wordmine <- wordmine[colnames(wordmine) != "digestive"]
wordmine <- t(wordmine)
rownames(wordmine)[grep("mouth", rownames(wordmine))] <- "oesophagus"
rownames(wordmine)[grep("biliary", rownames(wordmine))] <- "biliary_tract"


rownames(wordmine)[grep("hematopoietic", rownames(wordmine))] <- "haematopoietic_and_lymphoid_tissue"
wordmine["haematopoietic_and_lymphoid_tissue", ] <- wordmine["haematopoietic_and_lymphoid_tissue", ] + wordmine["blood", ]
wordmine <- wordmine[grep("blood", rownames(wordmine)) * -1, ]

wordmine["colorectal", ] <- wordmine["colorectal", ] + wordmine["colon", ] + wordmine["rectum", ]
wordmine <- wordmine[grep("colon", rownames(wordmine)) * -1, ]
wordmine <- wordmine[grep("rectum", rownames(wordmine)) * -1, ]

rownames(wordmine)[grep("gastrointestinal", rownames(wordmine))] <- "gastrointestinal_tract_(site_indeterminate)"
wordmine <- rbind(wordmine, wordmine["gastrointestinal_tract_(site_indeterminate)", ])
wordmine <- rbind(wordmine, wordmine["gastrointestinal_tract_(site_indeterminate)", ])
rownames(wordmine)[nrow(wordmine)] <- "small_intestine"
rownames(wordmine)[nrow(wordmine) - 1] <- "large_intestine"

wordmine <- melt(as.matrix(t(wordmine)))

save(wordmine, file = "./wordmine.RData")
