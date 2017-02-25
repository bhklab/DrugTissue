f <- file("./output/karyotype.txt")

if(!file.exists("./temp/figure7.xlsx"))
{
  source("./figure7p1.R")
}
edges <- read.xlsx("./temp/figure7.xlsx")

#edges$V2 <- gsub("_", " ", edges$V2)
colors <- c("blues-6-seq", "bugn-6-seq", "bupu-6-seq", "gnbu-6-seq", "greens-6-seq", "oranges-6-seq", "orrd-6-seq",
            "pubu-6-seq",
            "pubugn-6-seq",
            "purd-6-seq",
            "purples-6-seq",
            "rdpu-6-seq",
            "reds-6-seq",
            "ylgn-6-seq",
            "ylgnbu-6-seq",
            "ylorbr-6-seq",
            "ylorrd-6-seq")

c1 <- 1
c2 <- 1

boxes <- unique(c(edges$V1, edges$V2))
output <- c()
for(b in boxes)
{
  if(c2 > 6)
  {
    c2 <- 1
    c1 <- c1 + 1
  }
  if(c1 > length(colors)){c1 <- 1}
  
  cout <- paste(colors[c1], c2, sep = "-")
  output <- c(output, paste("chr -", b, b, 0, 100, cout, sep = " "))
  c2 <- c2 + 1
}

writeLines(output, f, useBytes = TRUE)
close(f)


output <- c()
f <- file("./output/edges.txt")
for(e in 1:nrow(edges))
{
  if(edges[e, "V3"] == 1)
  {
    n <- paste(edges[e, "V1"], 1, 25, edges[e, "V2"], 1, 25, "color=green_a3", sep = " ")
  }
  if(edges[e, "V3"] == 2)
  {
    n <- paste(edges[e, "V1"], 26, 50, edges[e, "V2"], 26, 50, "color=red_a3", sep = " ")
  }
  if(edges[e, "V3"] == 3)
  {
    n <- paste(edges[e, "V1"], 51, 75, edges[e, "V2"], 51, 75, "color=vdblue,z=80", sep = " ")
  }
  if(edges[e, "V3"] == 4)
  {
    n <- paste(edges[e, "V1"], 76, 100, edges[e, "V2"], 76, 100, "color=vlred_a3", sep = " ")
  }
  output <- c(n, output)
}
  
writeLines(output, f)
close(f)