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

cmb <- cmb <- cm[which(cm$value <= 0.05),]
cmb <- names(sort(-table(cmb$Var2))[c(1:3)])
maxdiff <- 0
maxdrug <- NULL

asd  <- table(cmw$Var2)
cms <- names(asd[asd == 26])

for(az in cms)
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
  pd <- na.omit(pd)
  if(median(pd$auc) > maxdiff)
  {
    maxdiff <- median(pd$auc)
    maxdrug <- az
  }
}

cmb <- c(cmb, "lapatinib", "ABT-888", maxdrug)

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
test[, "sig"] <- "no"

for(t in unique(test$tissue))
{
  if(length(grep(t, test$tissue)) < 10)
  {
    test <- test[which(test$tissue != t), ]
  }
}

for(drs in unique(test$drug))
{
  test2 <- test[which(test$drug == drs), ]

  for(qqq in unique(test2$tissue))
  {
    if(!is.na(combined1[qqq, drs]) && combined1[qqq, drs] <= 0.05)
    {
      test[which(test$drug == drs & test$tissue == qqq), "sig"] <- "yes"
    }
  }
}

test$drug <- as.factor(test$drug)

for(drgs in cmb)
{
  pdf(paste0("figure3 ",drgs, ".pdf"))
  test3 <- test[test$drug == drgs, ]
  print(
    ggplot(test3, aes(as.factor(test3$tissue), test3$auc)) + 
    geom_boxplot(outlier.colour = NA, aes(color = test3$sig)) + 
    labs(x = "tissue type", y = "recomputed AUC") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), 
          axis.text.y = element_text(size = 15), 
          axis.title=element_text(size=15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15)
          ) +
    scale_colour_manual(name = "significant", 
                        values = c("no" = "blue","yes" = "red")
                        ) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  )
  dev.off()
}

