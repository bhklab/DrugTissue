combinedcopy <- combined1
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
colnames(combinedcopy) <- toupper(gsub(badchars, "", colnames(combinedcopy)))
#combinedcopy <- combinedcopy[grep("breast_", rownames(combinedcopy)) * -1, ]

ccc <- combinedcopy

for(x in rownames(ccc))
{
  if(!(FALSE %in% is.na(ccc[x, ])))
  {
    ccc <- ccc[rownames(ccc) != x, ]
  }
}

for(x in colnames(ccc))
{
  if(!(FALSE %in% is.na(ccc[, x])))
  {
    ccc <- ccc[,colnames(ccc) != x ]
  }
}

cccm <- melt(as.matrix(ccc))
#cccm$value <- as.factor(cccm$value)
cccm <- na.omit(cccm)
#cccm$value <- as.factor(cccm$value)

cccm[, "drugc"] <- NA
for(a in 1:nrow(cccm))
{
  if(cccm[a, "Var2"] %in% toupper(gsub(badchars, "", breast)))
  {
    cccm[a, "drugc"] <- "breast"
  }
  if(cccm[a, "Var2"] %in% toupper(gsub(badchars, "", prostate)))
  {
    cccm[a, "drugc"] <- "prostate"
  }
  if(cccm[a, "Var2"] %in% toupper(gsub(badchars, "", colorectal)))
  {
    cccm[a, "drugc"] <- "colorectal"
  }
}

cccm <- na.omit(cccm)
cccm <- cccm[order(cccm$drugc), ]

jac <- combined1
jac <- log10(1/(jac + 1e-10))
jacc <- pvclust(t(jac), method.dist="cor", method.hclust="complete")


cccm$value <- log10(1/(cccm$value + 1e-10))

ggplot(cccm) + geom_tile(aes(fill = cccm$value, x = cccm$Var1, y = cccm$Var2)) + 
  scale_fill_continuous(name = "-log10 p value", low = "#0571b0", high = "#ca0020") + 
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 270, hjust = 1, size = 10), axis.text.y = element_text(size = 1)) +
  scale_x_discrete(limits = rownames(ccc)[jacc$hclust$order]) + 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))

#figure 6 supp

ccc <- combinedcopy

for(x in rownames(ccc))
{
  if(!(FALSE %in% is.na(ccc[x, ])))
  {
    ccc <- ccc[rownames(ccc) != x, ]
  }
}

for(x in colnames(ccc))
{
  if(!(FALSE %in% is.na(ccc[, x])))
  {
    ccc <- ccc[,colnames(ccc) != x ]
  }
}

cccm <- melt(as.matrix(ccc))
#cccm$value <- as.factor(cccm$value)
cccm <- na.omit(cccm)
#cccm$value <- as.factor(cccm$value)

cccm[, "drugc"] <- NA
for(a in 1:nrow(cccm))
{
  if(cccm[a, "Var2"] %in% toupper(gsub(badchars, "", breast)))
  {
    cccm[a, "drugc"] <- "breast"
  }
  if(cccm[a, "Var2"] %in% toupper(gsub(badchars, "", prostate)))
  {
    input <- data.frame(Var1 = cccm[a, "Var1"], Var2 = cccm[a, "Var2"], value = cccm[a, "value"], drugc = "prostate")
    cccm <- rbind(cccm, input)
  }
  if(cccm[a, "Var2"] %in% toupper(gsub(badchars, "", colorectal)))
  {
    input <- data.frame(Var1 = cccm[a, "Var1"], Var2 = cccm[a, "Var2"], value = cccm[a, "value"], drugc = "colorectal")
    cccm <- rbind(cccm, input)
  }
}

cccm <- na.omit(cccm)
cccm <- cccm[order(cccm$drugc), ]

jac <- combined1
jac <- log10(1/(jac + 1e-10))
jacc <- pvclust(t(jac), method.dist="cor", method.hclust="complete")


cccm$value <- log10(1/(cccm$value + 1e-10))

ggplot(cccm) + geom_tile(aes(fill = cccm$value, x = cccm$Var1, y = cccm$Var2)) + 
  scale_fill_continuous(name = "-log10 p value", low = "#0571b0", high = "#ca0020") + 
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 270, hjust = 1, size = 10), axis.text.y = element_text(size = 1)) +
  facet_grid(drugc ~ ., scales = "free_y") +
  scale_x_discrete(limits = rownames(ccc)[jacc$hclust$order]) + 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))
