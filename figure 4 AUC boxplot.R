options(stringsAsFactors = FALSE)

cm <- melt(as.matrix(combined1))
cm <- cm[cm$value != 2, ]
cm <- cm[order(cm$value), ]

# cmb <- cm[which(cm$value == 0),]
# cmb[, "size"] <- 0
# 
# for(x in 1:nrow(cmb))
# {
#   drg <- cmb[x, "Var2"]
#   for(y in mw)
#   {
#     if(cmb[x, "Var2"] %in% colnames(y) && cmb[x, "Var1"] %in% rownames(y))
#     {
#       
#       cmb[x, "size"] <- cmb[x, "size"] + y[cmb[x, "Var1"], which(colnames(y) == cmb[x, "Var2"])]
#     }
#   }
# }
# cmb <- cmb[order(-cmb$size),]
# cmb <- cmb[names(sort(-table(cma$Var2))[c(1:4)]), ]
# cmb[, "type"] <- "good"

cmb <- names(sort(-table(cma$Var2))[c(1:4)])

pd <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(pd) <- c("tissue", "drug", "auc")
for(az in cmb)
{
  for(z in c(CCLE, GDSC, CTRPv2, gCSI))
  {
    cle <- z@sensitivity$info[which(z@sensitivity$info$drugid == as.character(az)),]
    cle <- cle[, c("cellid", "drugid")]
    if(nrow(cle) > 0)
    {
      for(y in 1:nrow(cle))
      {
        cle[y, "tissue"] <- z@cell[cle[y, "cellid"], "tissueid"]
      }
      cle[, "auc"] <- z@sensitivity$profiles[rownames(cle), "auc_recomputed"]
      
      for(dp in 1:nrow(cle))
      {
        input <- data.frame(tissue = cle[dp, "tissue"], drug = az, auc = cle[dp, "auc"])
        pd <- rbind(pd, input)
      }
    }
  }
}

test <- pd
#test[, "label"] <- paste(test$drug, test$tissue , sep = " on \n ")

#test$label<- factor(test$label, unique(test$label))
#test$type <- as.factor(test$type)
test <- na.omit(test)
test <- test[test$tissue != "", ]
test[, "sig"] <- 0

for(t in unique(test$tissue))
{
  if(length(grep(t, test$tissue)) < 10)
  {
    test <- test[test$tissue != t, ]
  }
}

for(drs in unique(test$drug))
{
  test2 <- test[test$drug == drs, ]
  
  wilc <- matrix(nrow = max(table(test2$tissue)), ncol = length(unique(test2$tissue)))
  colnames(wilc) <- unique(test2$tissue)
  
  for(dd in unique(test2$tissue))
  {
    count <- 1
    for(ent in grep(dd, test2$tissue))
    {
      wilc[count, dd] <- as.numeric(test2[ent, "auc"])
      count <- count + 1
    }
  }
  op <- pairwise.wilcox.test(wilc, colnames(wilc), p.adjust.method = "none", alternative = "less")
  
  for(qqq in colnames(op$p.value))
  {
    qqqq <- op$p.value[, qqq]
    if(length(qqqq[which(qqqq <= 0.05)]) / length(na.omit(qqqq)) > 0)
    {
      test[test$drug == drs & test$tissue == qqq, "sig"] <- 1
    }
  }
}

test$sig <- factor(test$sig)

ggplot(test, aes(as.factor(test$tissue), test$auc)) + 
  geom_boxplot(outlier.colour = NA, aes(color = test$sig)) + 
  labs(x = "tissue type", y = "recomputed AUC") + 
  facet_grid(drug~.) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
  scale_colour_brewer(name = "significant", palette = "Set1", direction = -1)

