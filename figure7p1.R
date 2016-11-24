.libPaths("/rlibs/")


if(!file.exists("./temp/combined1.RData"))
{
  source("gsea_with_AUC.R")
} else 
{
  load("./temp/combined1.RData")
}

if(!file.exists("./temp/wordmine.RData"))
{
  source("wordmine.R")
} else 
{
  load("./temp/wordmine.RData")
}



options(stringsAsFactors = FALSE)
edges <- data.frame()
counter <- 1


# 1 = experimental
# 2 = clinical
# 3 = both 
# 4 = shows experimental clustering but lack of clinical importance 
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
combinedcopy <- combined1
colnames(combinedcopy) <- toupper(gsub(badchars, "", colnames(combinedcopy)))
#combinedcopy <- combinedcopy[grep("breast_", rownames(combinedcopy)) * -1, intersect(toupper(gsub(badchars, "", colorectal)), colnames(combinedcopy))]
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

for(x in 1:nrow(wordmine))
{
  if(wordmine$Var2[x] != "")
  {
    if(wordmine$value[x] >= 100)
    {
      y <- toupper(gsub(badchars, "", wordmine$Var1[x]))
      if(nrow(edges[edges$V1 == y & edges$V2 == wordmine$Var2[x],]) == 0)
      {
        if(y %in% colnames(combinedcopy))
        {
          if(combinedcopy[wordmine$Var2[x],y] <= 0.05 && length(combinedcopy[wordmine$Var2[x], y]) > 0 && !is.na(combinedcopy[wordmine$Var2[x], y]))
          {
            edges[counter, 1] <- y
            edges[counter, 2] <- as.character(wordmine[x, "Var2"])
            message(wordmine[x, "Var2"])
            edges[counter, 3] <- 3
            counter <- counter + 1
          }
          else
          {
            edges[counter, 1] <- y
            edges[counter, 2] <-as.character(wordmine[x, "Var2"])
            message(wordmine[x, "Var2"])
            edges[counter, 3] <- 2
            counter <- counter + 1
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
    if(grepl("breast", a))
    {
      a <- "breast"
    }
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

edgescopy <- edges

edges <- edgescopy

for(a in 1:nrow(edges))
{
  if(edges[a, 3] == 1)
  {
    nsamples <- 0
    nmatch <- 0
    for(b in c(gCSI, CCLE, GDSC1000))
    {
      for(c in grep(edges[a, 1], toupper(gsub(badchars, "", b@sensitivity$info$drugid))))
      {
        ttype <- b@cell[b@sensitivity$info[c, "cellid"], "tissueid"]
        if(!is.na(ttype) && ttype == edges[a,2] && length(ttype) != 0)
        {
          auc <- b@sensitivity$profiles[rownames(b@sensitivity$info)[c], "auc_recomputed"]
          if(length(auc) != 0 && !is.na(auc))
          {
            nsamples <- nsamples + 1
            if(auc >= 0.2)
            {
              nmatch <- nmatch + 1
            }
          }
        }
      }
    }
    b <- CTRPv2
    {
      for(c in grep(edges[a, 1], toupper(gsub(badchars, "", b@sensitivity$info$drugid))))
      {
        ttype <- b@cell[b@sensitivity$info[c, "cellid"], "tissueid"]
        if(!is.na(ttype) && ttype == edges[a,2] && length(ttype) != 0)
        {
          auc <- b@sensitivity$profiles[rownames(b@sensitivity$info)[c], "auc_recomputed"]
          if(length(auc) != 0 && !is.na(auc))
          {
            nsamples <- nsamples + 1
            if(auc >= 0.2)
            {
              nmatch <- nmatch + 1
            }
          }
        }
      }
    }
    if(nsamples == 0)
    {
      edges[a, 3] <- NA
    }
    else if(nmatch / nsamples < 0.1)
    {
      edges[a,3] <- 4
    }
  }
}

edges <- na.omit(edges)
write.xlsx(edges, file = "./temp/figure7.xlsx")
