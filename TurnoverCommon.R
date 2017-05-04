DataDir = "~/Google Drive/CoChIPAnalysis/Turnover"
#setwd("~/Dropbox/coChIP-scripts/")
library(abind)
source("CoChip-Functions.R")
source("NucAtlas.R")
source("Multiplot.R")

## prepere SGD/nuc structures
GetSGDandNucs()
source("Main-Functions.R")

if( !exists("Data") ) {
  print("Aggregating Turnover Data")
  
  NucsFN = paste0(DataDir, "/Nucs.rdata")
  Nucs <- readRDS(NucsFN)
  colnames(Nucs) = 1:(dim(Nucs)[2])
  Names = rownames(Nucs)
  Strain <- factor(unlist(lapply(strsplit(Names,"-"), function(x) paste0(x[[1]],"-",x[[2]]))))
  Ab1 <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[3]])))
  Ab2 <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[4]])))
  
  Time.0CC = 125
  Time.1CC = 150
  Time.2CC = 175
  
  NucsTime = list()
  Time = sub("Raf-Gal", "", levels(Strain))
  Time = sub("Gal-Glu0", Time.0CC, Time)
  Time = sub("Gal-Glu1CC", Time.1CC, Time)
  Time = sub("Gal-Glu2CC", Time.2CC, Time)
  names(Time) = levels(Strain)
  
  for( s in c("Flag", "Myc")) {
    NucsTime[[s]] = list()
    for( a in c("Input", "H3K4me3", "H3K36me3", "H3K79me3", "Beads")) {
      TargetSet = Ab2 == a & Ab1 == s
      Val = Nucs[TargetSet,]
      rownames(Val) = Time
      NucsTime[[s]][[a]] = Val
    }
  }
  
  Data = abind(lapply(NucsTime, function(x) abind(x,along=-1)), along=-1)
}  
