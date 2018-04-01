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


TargetMods = c("H3K4me3")
mod.nucs = ReadHistoneModeAtlas("Weiner-HisMod.csv", NucRegions, TargetMods )
m = mcols(mod.nucs)
weiner.nucs = do.call(rbind,lapply(TargetMods, function(mod) as.vector(exp(as.matrix(m["Input"]) + as.matrix(m[mod])))))
weiner.nucs = rbind(weiner.nucs, as.vector(exp(as.matrix(m["Input"]))))
rownames(weiner.nucs) = c(TargetMods, "Input")
occ.weiner = pmin(weiner.nucs["Input",] / quantile(weiner.nucs["Input",], 0.99, na.rm = TRUE),1)

occ.CC0 = Data["Flag", "Input", "125",]
occ.CC0 = pmin(occ.CC0 / quantile(occ.CC0, 0.99, na.rm = TRUE), 1)

occ.M10 = Data["Myc", "Input", "10",]
occ.M10 = pmin(occ.M10 / quantile(occ.M10, 0.99, na.rm = TRUE), 1)

occ.emp = .5*(occ.CC0+occ.M10)

p.Occ = "Occ"
p.Rate = "Rate"
p.Err = "Err"
p.T0 = "T0"
p.Flag = "Flag"
p.Myc = "Myc"

p.Names = c(p.Occ, p.Rate, p.Err, p.T0, p.Flag, p.Myc)

ModelParamaterVector <- function() {
  p = rep(1,length(p.Names))
  names(p) = p.Names
  p
}

ReadsByTime = function( t, p, w = rep(1, length(t)), Flag = TRUE ) {
  c = p[p.Occ]
  if( Flag )  {
    c = c * p[p.Flag]
  } else
    c = c * p[p.Myc]
  r = p[p.Rate]
  e = p[p.Err]
  t0 = p[p.T0]
  T = HingeApprox(t-t0)
  if( Flag ) {
    w*(c*(1-exp(-r*T)) + e)
  } else {
    w*(c*exp(-r*T) + e)
  }
}

ReadsByTime.Grad = function( t,p, w = rep(1, length(t)), Flag = TRUE ) {
  c = p[p.Occ]
  if( Flag )  {
    c = c * p[p.Flag]
  } else
    c = c * p[p.Myc]
  r = p[p.Rate]
  e = p[p.Err]
  t0 = p[p.T0]
  T = HingeApprox(t-t0)
  gr = matrix( 0, nr = length(p), nc = length(t))
  rownames(gr) = names(p)
  colnames(gr) = t
  if( Flag ) {
    gr[p.Occ,] = w*(1-exp(-r*T))*c/p[p.Occ]
    gr[p.Flag,] = w*(1-exp(-r*T))*p[p.Occ]
    gr[p.Rate,] = w*c*T*exp(-r*T)
    gr[p.T0,] = -w*c*r*HingeApprox.Grad(t-t0)*exp(-r*T)
  } else {
    gr[p.Myc,] = -w*(exp(-r*T))*p[p.Occ]
    gr[p.Occ,] = w*(exp(-r*T))*c/p[p.Occ]
    gr[p.Rate,] = -w*c*T*exp(-r*T)
    gr[p.T0,] = w*c*r*HingeApprox.Grad(t-t0)*exp(-r*T)
  }
  gr[p.Err,] = w
  gr
}


HingeApprox = function(x,alpha=10){
  x+log(1+exp(-alpha*x))/alpha
}

HingeApprox.Grad = function(x,alpha=10){
  1-alpha*exp(-alpha*x)/(alpha*(1+exp(-alpha*x)))
}

FitNucsFN = paste0(DataDir, "/FitNucs.rdata")
Thresholds = sapply(1:dim(Data)[3], function(i) quantile(Data["Flag", "Input",i,], .995))
NucIDs = which(apply(sapply(1:dim(Data)[3], function(i) Data["Flag", "Input",i,] <= Thresholds[i] ), 1, all))



if( file.exists(FitNucsFN)) {
  FitNucs = readRDS(FitNucsFN)
  
  params.Occ = sapply(FitNucs,function(x) x$fix.param[p.Occ])
  params.Rate = sapply(FitNucs,function(x) x$fix.param[p.Rate])
}

