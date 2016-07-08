library(reshape2)

wordmine <- melt(as.matrix(wordmine))
wordmine <- na.omit(wordmine)
wordmine$value <- log10(wordmine$value)
hist(wordmine$value, xlab = "log10 adjusted hits", main = "")

plot(hclust(dist(t(ccc3))))
pvclust::pvclust(combined1, method.hclust = "complete", method.dist = "correlation")
plot(hclust(vegdist(combined1, method = "jaccard")))
