options(stringsAsFactors = FALSE)

.libPaths("/rlibs/")
library(PharmacoGx)

if(!file.exists("./temp/combined1.RData"))
{
  source("./gsea_with_AUC.R")
} else {
  load("./temp/combined1.RData")
}

CCLE <- downloadPSet("CCLE_2013")
GDSC100 <- downloadPSet("GDSC1000")
gCSI <- downloadPSet("gCSI")
CTRPv2 <- downloadPSet("CTRPv2")


f3a <- data.frame(matrix(ncol = 1, nrow = nrow(combined1)))
rownames(f3a) <- rownames(combined1)
colnames(f3a) <- "value"
for(x in rownames(combined1))
{
  f3a[x, 1] <- length(combined1[!is.na(combined1) & combined1 <= 0.05 & rownames(combined1) == x])
}

order <- rownames(f3a)[order(-f3a$value)]

pdf("figure4a.pdf")
ggplot(f3a, aes(x = rownames(f3a), y = f3a$value)) + 
  geom_bar(stat = "identity", position="dodge") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.title=element_text(size=10)) +
  guides(fill=guide_legend(title="")) +
  scale_x_discrete(limits = order) + 
  labs(x = NULL, y = "number of significant interactions") +
  theme(legend.text=element_text(size=10)) +
  coord_flip()
dev.off()


f3b <- data.frame(matrix(ncol = 2, nrow = ncol(combined1)))
rownames(f3b) <- as.factor(colnames(combined1))
colnames(f3b) <- c("drug", "value")
f3b$drug <- colnames(combined1)

for(x in colnames(combined1))
{
  xd <- na.omit(combined1[,x])
  f3b[x, "value"] <- length(xd[xd <= 0.05])
}


orv <- order(-f3b$value)
top3 <- rownames(f3b)[orv[1:4]]
bot3 <- rownames(f3b)[grep(0, f3b$value)]

orv <- orv[1:nrow(combined1)]
f3b <- f3b[orv, ]

pdf("figure4b.pdf")
ggplot(f3b, aes(x = rownames(f3b), y = f3b$value)) + 
  geom_bar(stat = "identity", position="dodge") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.title=element_text(size=10)) +
  guides(fill=guide_legend(title="")) +
  scale_x_discrete(limits = rownames(f3b)) + 
  labs(x = NULL, y = "number of significant interactions") +
  theme(legend.text=element_text(size=10)) +
  coord_flip()
dev.off()


cmb <- unique(c(top3, bot3))

pd <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(pd) <- c("tissue", "drug", "auc")

for(az in cmb)
{
  for(z in c(CCLE, GDSC1000, CTRPv2, gCSI))
  {
    cle <- z@sensitivity$info[which(z@sensitivity$info$drugid == az),]
    cle <- cle[, c("cellid", "drugid")]
    if(nrow(cle) > 0)
    {
      for(y in 1:nrow(cle))
      {
        if(z@cell[cle[y, "cellid"], "tissueid"] == "esophagus" && !is.na(z@cell[cle[y, "cellid"], "tissueid"]))
        {
          cle[y, "tissue"] <- "oesophagus"
        }else
        {
          cle[y, "tissue"] <- z@cell[cle[y, "cellid"], "tissueid"]
        }
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
  pdf(paste0("figure3 ",gsub(badchars, "", drgs), ".pdf"))
  test3 <- test[test$drug == drgs, ]
  print(
    ggplot(test3, aes(as.factor(test3$tissue), test3$auc)) + 
      geom_boxplot(outlier.colour = NA, aes(color = test3$sig)) + 
      labs(x = "tissue type", y = "recomputed AUC") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
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



