linenums <- as.data.frame(matrix(ncol = 4, nrow = nrow(combined1)))
rownames(linenums) <- rownames(combined1)
colnames(linenums) <- c("CCLE", "GDSC", "CTRPv2", "gCSI")
linenums[is.na(linenums)] <- 0

for(y in 1:nrow(CCLE@cell))
{
  if(length(grep(CCLE@cell$cellid[y], rownames(na.omit(CCLE@sensitivity$profiles)))) > 0)
  {
    if(!is.na(CCLE@cell$tissueid[y]))
    {
      linenums[CCLE@cell$tissueid[y], "CCLE"] <- linenums[CCLE@cell$tissueid[y], "CCLE"] + 1
    }
  }
}
for(y in 1:nrow(GDSC@cell))
{
  if(length(grep(GDSC@cell$cellid[y], rownames(na.omit(GDSC@sensitivity$profiles)))) > 0)
  {
    if(!is.na(GDSC@cell$tissueid[y]))
    {
      linenums[GDSC@cell$tissueid[y], "GDSC"] <- linenums[GDSC@cell$tissueid[y], "GDSC"] + 1
    }
  }
}
for(y in 1:nrow(CTRPv2@cell))
{
  if(length(grep(CTRPv2@cell$cellid[y], rownames(na.omit(CTRPv2@sensitivity$profiles)))) > 0)
  {
    if(!is.na(CTRPv2@cell$tissueid[y]))
    {
      linenums[CTRPv2@cell$tissueid[y], "CTRPv2"] <- linenums[CTRPv2@cell$tissueid[y], "CTRPv2"] + 1
    }
  }
}
for(y in 1:nrow(gCSI@cell))
{
  if(length(grep(gCSI@cell$unique.id[y], rownames(na.omit(gCSI@sensitivity$profiles)))) > 0)
  {
    if(!is.na(gCSI@cell$tissueid[y]))
    {
      linenums[gCSI@cell$tissueid[y], "gCSI"] <- linenums[gCSI@cell$tissueid[y], "gCSI"] + 1
    }
  }
}


linenums <- linenums[grep("breast_", rownames(linenums)) * -1, ]
linenums <- na.omit(linenums)

for(x in 1:nrow(linenums))
{
  linenums[x, "median"] <- median(as.numeric(linenums[x, 1:4]))
}

linenums <- linenums[order(linenums$median), ]

linenums <- linenums[, -5]
order <- rownames(linenums)

linenums <- melt(as.matrix(linenums))

colnames(linenums) <- c("tissue type", "dataset", "value")

linenums$`tissue type` <- gsub("_", " ", linenums$`tissue type`)

linenums$`tissue type` <- as.factor(linenums$`tissue type`)

pdf("figure2a.pdf")
ggplot(linenums, aes(x = linenums$`tissue type`, y = linenums$value, fill=linenums$dataset)) + 
  geom_bar(stat = "identity", position="dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.title=element_text(size=10)) +
  guides(fill=guide_legend(title="")) +
  scale_x_discrete(limits = gsub("_", " ", order)) + 
  labs(x = NULL, y = "number of cell lines") +
  theme(legend.text=element_text(size=10))
dev.off()
  