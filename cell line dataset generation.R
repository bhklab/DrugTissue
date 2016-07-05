linenums <- linenums[grep("breast_", rownames(linenums)) * -1, ]

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

ggplot(linenums, aes(x = linenums$`tissue type`, y = linenums$value, fill=linenums$dataset)) + 
  geom_bar(stat = "identity", position="dodge") + 
  theme(axis.text.x = element_text(angle = 270, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.title=element_text(size=10)) +
  guides(fill=guide_legend(title="")) +
  scale_x_discrete(limits = gsub("_", " ", order)) + 
  labs(x = NULL, y = "number of cell lines") +
  theme(legend.text=element_text(size=10))
  