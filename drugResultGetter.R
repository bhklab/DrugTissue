getDrugTests <- function(drug)
{
  drug <- drug
  url <- paste("https://clinicaltrials.gov/search?term=", drug, sep="")
  url <- paste(url, "&studyxml=true", sep="")
  
  cellannotation <- read.csv(file.path("PharmacoGx-Private/inst/extdata", "cell_annotation_all.csv"), sep=",", comment.char="#") 
  tissues <- cellannotation[, grep("tissue", colnames(cellannotation))]
  to <- vector()
  
  for(i in 1:ncol(tissues))
  {
    to <- c(to, tissues[, i])
  }
  
  to <- c(to, "colorectal")
  
  to <- to[to != ""]
  to <- tolower(to)
  to <- unique(to)
  to <- gsub(" ", "", to)
  to <- gsub("_", " ", to)
  
  for(x in 1:length(to))
  {
    if(nchar(to[x]) < 4)
    {
      to[x] <- NA
    }
    else
    {
      if(length(grep(" ", to[x])) > 0)
      {
        to[x] <- unlist(strsplit(to[x], split=" "))
      }
    }
  }
  to <- na.omit(to)
  if(!file.exists(paste(drug, ".zip", sep="")))
  {
    download.file(url=url, destfile = paste(drug, ".zip", sep=""), mode="wb", quiet = TRUE)
  }
  dir.create(drug, showWarnings=FALSE, recursive=TRUE) 
  if(file.exists(paste(drug, ".zip", sep="")))
  {
    unzip(zipfile=paste(drug, ".zip", sep=""), exdir = drug, overwrite = FALSE)
  }
  else
  {
    return(NA)
  }
  
  
  xmls <- list.files(drug)
  trialnums <- length(xmls)
  
  output <- vector()
  zz <- 1;
  for(x in xmls)
  {
      t <- XML::xmlParse(file.path(drug,x))
      u <- XML::xmlSApply(XML::xmlRoot(t), function(x) XML::xmlSApply(x, XML::xmlValue))
      v <- data.frame(t(u),row.names=NULL)
      
      
      a <- v[, "brief_summary"]
      a <- gsub("\r?\n|\r", " ", a)
      a <- unlist(strsplit(a, split=" "))
      b <- grep("cancer", a, ignore.case = TRUE, fixed = TRUE)
      for(i in b)
      {
        if(a[i-1] %in% to)
        {
          #print(a[i-1])
          output[zz] <- gsub("[[:punct:]]", "", a[i-1])
          zz <- zz + 1
        }
      }
      
      for(y in grep("keyword", colnames(v), ignore.case = TRUE))
      {
        if(length(grep("cancer", v[, y]$keyword[[1]], ignore.case = TRUE, fixed = TRUE)) > 0)
        {
          
          l<- gsub("[[:punct:]]", "", v[, y]$keyword[[1]])
          l <- unlist(strsplit(a, split=" "))
          m <- grep("cancer", l, ignore.case = TRUE, fixed = TRUE)
          for(i in m)
          {
            if(l[i-1] %in% to)
            {
              #print(l[i-1])
              output[zz] <- gsub("[[:punct:]]", "", l[i-1])
              zz <- zz + 1
            }
          }
        }
      }
      
      for(y in grep("condition", colnames(v), ignore.case = TRUE))
      {
        if(length(grep("cancer", v[, y]$keyword[[1]], ignore.case = TRUE, fixed = TRUE)) > 0)
        {
          
          l<- gsub("[[:punct:]]", "", v[, y]$keyword[[1]])
          l <- unlist(strsplit(a, split=" "))
          m <- grep("cancer", l, ignore.case = TRUE, fixed = TRUE)
          for(i in m)
          {
            if(l[i-1] %in% to)
            {
              #print(l[i-1])
              output[zz] <- gsub("[[:punct:]]", "", l[i-1])
              zz <- zz + 1
            }
          }
        }
      }
      #output <- unique(output)
  }
  th <- output
  
  output <- th
  output <- gsub("cancers", "cancer", output)
  output <- tolower(output)
  output <- na.omit(output)
  
  
  if(length(output) > 100)
  {
    for(i in unique(output))
    {
      if(length(grep(i, output)) < 2)
      {
        for(j in 1:length(output))
        {
          if(output[j] == i)
          {
            output[j] <- NA
          }
        }
        output <- na.omit(output)
      }
    }
  }
  
  message(drug)
  ooo <- data.frame(matrix(ncol = 2, nrow = length(unique(output))))
  ooo[, 2] <- unique(output)
  if(nrow(ooo) > 0)
  {
    for(i in 1:nrow(ooo))
    {
      if(!is.na(ooo[i,2]))
      {
        ooo[i,1] <- length(grep(ooo[i, 2], output))
      }
    }
  }
  
  if(nrow(ooo) == 0)
  {
    return(NA)
  }
  else
  {
    return(ooo)
  }
  
}
  