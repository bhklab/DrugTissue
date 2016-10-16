boxp <- data.frame(matrix(nrow = nrow(combined1), ncol = 2))
colnames(boxp) <- c("percentage", "size")
rownames(boxp) <- rownames(combined1)

linenumscp <- melt(as.matrix(linenums))

for(c in rownames(combined1))
{
  boxp[c, "percentage"] <- length(which(combined1[c,] <= 0.05)) / length(na.omit(t(combined1[c,])))
  boxp[c, "size"] <- sum(linenumscp[grep(c, linenumscp$Var1, fixed = TRUE), "value"])
}

boxp[, "percentage"] <- as.numeric(boxp[, "percentage"] * 100)
boxlevel <- rownames(boxp)[order(boxp[,"percentage"])]

pdf("figure3a.pdf")
ggplot(data = boxp, aes(x = rownames(boxp), y = percentage)) + 
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
  scale_x_discrete(limits=boxlevel) +
  labs(x = "tissue type")
dev.off()

#spearman correlation
boxp <- boxp[grep("breast_", rownames(boxp)) * -1, ]
spear <- cor.test(x = boxp$percentage, y = boxp$size, method = "spearman")

