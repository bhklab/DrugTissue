library(reshape2)

wordmine <- melt(as.matrix(wordmine))
wordmine <- na.omit(wordmine)
wordmine$value <- log10(wordmine$value)
hist(wordmine$value, xlab = "log10 adjusted hits", main = "")
