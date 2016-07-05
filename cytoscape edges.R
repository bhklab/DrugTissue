edges <- data.frame()
counter <- 1

# 1 = experimental
# 2 = clinical
# 3 = both 
# 4 = wordmined

combinedcopy <- combined1
combinedcopy[is.na(combinedcopy)] <- 0
wordmine[is.na(wordmine)] <- 0
colnames(combinedcopy) <- toupper(gsub(badchars, "", colnames(combinedcopy)))
combinedcopy <- combinedcopy[grep("breast_", rownames(combinedcopy)) * -1, ]

for(a in unique(breast))
{
    a <- toupper(gsub(badchars, "", a))
    if(a %in% colnames(combinedcopy))
    {
      if(combinedcopy["breast", a] <= 0.05 && !is.na(combinedcopy["breast",a]))
      {
        edges[counter, 1] <- a
        edges[counter, 2] <- "breast"
        edges[counter, 3] <- 3
        counter <- counter + 1
      }
      else
      {
        edges[counter, 1] <- a
        edges[counter, 2] <- "breast"
        edges[counter, 3] <- 2
        counter <- counter + 1
      }
    }
}

for(a in unique(colorectal))
{
  a <- toupper(gsub(badchars, "", a))
  if(a %in% colnames(combinedcopy))
  {
    if(combinedcopy["large_intestine", a] <= 0.05 && !is.na(combinedcopy["large_intestine",a]))
    {
      edges[counter, 1] <- a
      edges[counter, 2] <- "large_intestine"
      edges[counter, 3] <- 3
      counter <- counter + 1
    }
    else
    {
      edges[counter, 1] <- a
      edges[counter, 2] <- "large_intestine"
      edges[counter, 3] <- 2
      counter <- counter + 1
    }
  }
}

for(a in unique(prostate))
{
  a <- toupper(gsub(badchars, "", a))
  if(a %in% colnames(combinedcopy))
  {
    if(combinedcopy["prostate", a] <= 0.05 && !is.na(combinedcopy["prostate",a]))
    {
      edges[counter, 1] <- a
      edges[counter, 2] <- "prostate"
      edges[counter, 3] <- 3
      counter <- counter + 1
    }
    else
    {
      edges[counter, 1] <- a
      edges[counter, 2] <- "prostate"
      edges[counter, 3] <- 2
      counter <- counter + 1
    }
  }
}

for(x in rownames(wordmine))
{
  for(y in colnames(wordmine))
  {
    if(x != "")
    {
      if(wordmine[x,y] >= 10 && x != "")
      {
        y <- toupper(gsub(badchars, "", y))
        if(nrow(edges[edges$V1 == y & edges$V2 == x,]) == 0)
        {
          if(y %in% colnames(combinedcopy))
          {
            if(combinedcopy[x, y] <= 0.05 && length(combinedcopy[x, y]) > 0 && !is.na(combinedcopy[x, y]))
            {
              edges[counter, 1] <- y
              edges[counter, 2] <- x
              edges[counter, 3] <- 3
              counter <- counter + 1
            }
            else
            {
              edges[counter, 1] <- y
              edges[counter, 2] <- x
              edges[counter, 3] <- 2
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
}

edgescopy <- edges
edges <- edgescopy

for(a in rownames(combinedcopy))
{
  for(b in colnames(combinedcopy))
  {
    if(combinedcopy[a,b] <= 0.05 && !is.na(combinedcopy[a,b]))
    {
      if(b %in% edges$V1 && nrow(edges[edges$V1 == b & edges$V2 == a,]) == 0)
      {
        edges[counter, 1] <- b
        edges[counter, 2] <- a
        edges[counter, 3] <- 1
        counter <- counter + 1
      }
    }
  }
}
