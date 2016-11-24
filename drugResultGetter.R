getDrugTests <- function(drug, phase3 = FALSE)
{
  .libPaths("/rlibs/")
  if(!require(XML))
  {
    install.packages("XML")
  }
  
  library(XML)
  options(stringsAsFactors = FALSE)
  drug <- drug
  path <- paste("./drugdata", drug, "", sep = "/")
  if(phase3)
  {
    url <- paste("https://clinicaltrials.gov/search?term=", drug,"+AND+phase+3", sep="")
  }else
  {
    url <- paste("https://clinicaltrials.gov/search?term=", drug,sep="")
  }
  url <- paste(url, "&studyxml=true", sep="")
  
  cellannotation <- read.csv(file.path("~/Documents", "cell_annotation_all_new.csv"), sep=",", comment.char="#") 
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
  to <- na.omit(to)
  for(x in 1:length(to))
  {
    if(nchar(to[x]) < 4 || is.na(to[x]))
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
  
  if(!file.exists("drugdata"))
  {
    dir.create("drugdata", showWarnings=FALSE, recursive=TRUE) 
  }
  
  if(!file.exists(path))
  {
    dir.create(path, showWarnings = FALSE, recursive= TRUE)
  }
 
  if(!file.exists(paste(path, drug, ".zip", sep="")))
  {
    download.file(url=url, destfile = paste(path, drug, ".zip", sep=""), mode="wb", quiet = TRUE)
  }
  
  if(file.exists(paste0(path, drug, ".zip", sep="")) && file.info(paste0(path, drug, ".zip", sep=""))$size > 0)
  {
    unzip(zipfile=paste(path, drug, ".zip", sep=""), exdir = path, overwrite = FALSE)
  }else
  {
    return(NA)
  }
  
  xmls <- list.files(path, pattern = "\\.xml$")
  trialnums <- length(xmls)
  
  output <- vector()
  zz <- 1;
  for(x in xmls)
  {
      t <- XML::xmlParse(file.path(path,x))
      u <- XML::xmlSApply(XML::xmlRoot(t), function(x) XML::xmlSApply(x, XML::xmlValue))
      v <- data.frame(t(u),row.names=NULL)
      
      if((phase3 && 3 %in% v$phase$phase[[1]]) || !phase3)
      {
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
  