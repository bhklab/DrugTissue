#figure 5
c3 <- combined1[is.na(combined1)] <- 2
jac.dist <- vegdist(combined1, method = "jaccard")
plot(hclust(jac.dist))


#figure 5.1
jac <- combined1
jacc <- pvclust(t(jac), method.dist="cor", method.hclust="complete")
plot(jacc)

#figure 5.2
for(x in 1:nrow(jac))
{
  for(y in 1:ncol(jac))
  {
    if(!is.na(jac[x,y]))
    {
      if(jac[x,y] > 0.05)
      {
        jac[x,y] <- 0
      }
      else
      {
        jac[x,y] <- 1
      }
    }
  }
}

y <- 1
while(y <= ncol(jac))
{
  if(NA %in% jac[,y] || length(grep(0, jac[,y])) == length(jac[,y]))
  {
    jac <- jac[, y * -1]
    y <- y - 1
  }
  y <- y + 1
}

for(x in rownames(jac))
{
  if(length(grep(0, jac[x,])) == ncol(jac))
  {
    jac <- jac[rownames(jac) != x, ]
  }
}

jac.dist <- vegdist(jac, method = "jaccard")
plot(hclust(jac.dist))

