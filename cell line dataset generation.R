for(x in 1:nrow(linenums))
{
  linenums[x, "median"] <- median(as.numeric(linenums[x, 1:4]))
}

linenums <- linenums[order(linenums$median), ]

linenums <- linenums[, -5]
order <- rownames(linenums)

linenums <- melt(as.matrix(linenums))

colnames(linenums) <- c("tissue type", "dataset", "value")

ggplot(linenums, aes(x = linenums$`tissue type`, y = linenums$value, fill=linenums$dataset)) + 
  geom_bar(stat = "identity", position="dodge") + 
  theme(axis.text.x = element_text(angle = 270, hjust = 1)) +
  scale_x_discrete(limits = order) +
  labs(x = "tissue types", y = "number of cell lines")