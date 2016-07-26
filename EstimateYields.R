WorkDir = "~/Dropbox/CoChIP/coChIP-Analysis"
DataDir = getwd()
OrigDir = getwd()
BamDir = paste0(DataDir,"/BAM")
setwd(OrigDir)

#DataDir =  "~/Data/CoChIP/160124/RPD3/"
extendDir <- function(x) {
  a = substr(x,1,1)
  if( a == "/" || a == "~")
    return(x)
  else
    return(paste0(OrigDir,"/",x))
}

#K4 experiment
LoadingRatio = list( H4 = 6.4,
                     K4me1 = 7.2,
                     K4me2 = 3.8,
                     K4me3 = 3.5,
                     K4ac = 17,
                     K4un = 26,
                     Input = 1)

#RPD3 experiment
LoadingRatio = list( H4 = 7.9, 
                     K4me3 = 4.9, 
                     K36me3 = 6.6, 
                     K79me3 = 2.1, 
                     K18ac = 24.8, 
                     K56ac = 3.7, 
                     Input = 1
  )

print("Initializing")
setwd(WorkDir)
suppressMessages(source("CoChip-Functions.R"))
suppressMessages(source("NucAtlas.R"))
source("MixPoisson.R")
GetSGDandNucs()


suppressMessages(library("optparse"))

option_list = list(
  make_option("--maxfraglen", type = "integer", default = 600,
              help="maximal fragment length"),
  make_option("--minfraglen", type = "integer", default = 50,
              help="minimum fragment length"),
  make_option(c("-d", "--datadir"), type = "character", default = NULL,
              help="data directory"),
  make_option(c("-b", "--bamdir"), type="character", default=NULL,
              help="location of BAM files" )
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments = TRUE);
## prepere SGD/nuc structures


params <- ccParams()
#params$MaxFragLen = 220
#params$MinFragLen = 50
params$MaxFragLen = opt$options$maxfraglen
params$MinFragLen = opt$options$minfraglen
params$cov = FALSE
params$cen = FALSE
params$meta = FALSE
params$nuc = FALSE

combine.tables <- function( t1, t2 ) {
  V = unique(c(levels(t1[,1]), levels(t2[,1])))
  t = data.frame( V1 = V, Freq = rep(0,length(V)))
  for( i in 1:length(V) )
    t[i,2] = sum(c(t1[t1[,1] == V[[i]],2], t2[t2[,1] == V[[i]],2]))
  
  t
}

ComputeTopRatio <- function(m1, m2, threshold = .25) {
  t1 = quantile(m1, threshold)
  t2 = quantile(m2, threshold)
  T1 = quantile(m1, .99)
  T2 = quantile(m2, .99)
  B = (m1 >= t1) & (m2 >= t2) & (m1 < T1) & (m2 < T2)
  #  print(length(which(B)))
  if( any(B) ) {
    data = data.frame(x = m1[B], y = m2[B])
    #    lm = lm(y~x-1, data = data)
    tryCatch({
      rlm = rlm(y~x-1, data = data)
      return(coef(rlm))
    }, error = function(e) print(e) )
    
    #    print(coef(rlm))
  } 
  return(NA)
}

if( !is.null(opt$options$datadir) )
  DataDir = extendDir(opt$options$datadir)
params$DataDir = DataDir

#Files = c("bar1-K18ac-Input", "bar1-K4me3-Input","bar1-K18ac-K18ac", "bar1-K4me3-K4me3")

Files = list.files( params$DataDir, "rdata$", full.names = F )
Names = Files
Names = sub("bar1", "WT", Names)
Names = sub(".rdata", "", Names)
Names = sub("_merge", "", Names)

Strain <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[1]])))
Ab1 <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[2]])))
Ab2 <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[3]])))

DoProcessFile <- function( x ) {
  f = BaseFileName(x)
  d = ccProcessFile(paste0(BamDir,"/",x), param = params)
  dat = d$GR
  dat = dat[width(dat) <= params$MaxFragLen & width(dat) > params$MinFragLen]
  dat.uniq = unique(dat)
  # count duplicate segments
  dups <- countOverlaps(dat.uniq, dat, type = "equal")
  dups = as.data.frame(table(dups))
  list(dups = dups, Nuc = d$nuc)
}

## do the work...
Dat = lapply(Files, DoProcessFile)
Hists = lapply(Dat, function(x) x$dups)
N.uniq = sapply(Hists, function(x) sum(x$Freq) )

Nucs = do.call(rbind, lapply(Dat, function(x) x$Nuc))
rownames(Nucs) = Names

CombinedHist = list()
for( i in 1:length(Files)) {
  ab2 = as.character(Ab2[[i]])
  if( ab2 %in% names(CombinedHist)) {
    CombinedHist[[ab2]] = combine.tables(CombinedHist[[ab2]], Hists[[i]])
  } else
    CombinedHist[[ab2]] = Hists[[i]]
}


Summary = data.frame(Poisson.Lambda = c(), Poisson.total = c(), Poisson.yield = c(), Mix.Lambda = c(), Mix.Total = c(), Mix.yield = c() )

for( a in names(CombinedHist))
{
  dupStat = CombinedHist[[a]]$Freq
  names(dupStat) = CombinedHist[[a]]$V1
  Nuniq = sum(dupStat)
  lambda1 = FitNonZeroPoisson(dupStat)
  Total1 = Nuniq /(1 - exp(-lambda1))
  Summary[a,"Poisson.Lambda"] = lambda1
  Summary[a,"Poisson.total"] = Total1
  Summary[a,"Poisson.yield"] = Nuniq/Total1
  
  p = FitMixPoisson(dupStat)
  lambda = p$lambda * sum(unlist(lapply(1:length(p$p), function (i) i*p$p[[i]])))
  Total = Nuniq /(1 - MixPoisson(0, p$lambda, p$p))
  Summary[a,"Mix.Lambda"] = lambda
  Summary[a,"Mix.Total"] = Total
  Summary[a,"Mix.yield"] = Nuniq/Total
}

print(Summary)

AbYield = list( Input = 1)
for( a in levels(Ab1) )
{
  yield = list()
  ryield = list()
  for( s in levels(Strain) )
  {
    Input = which(Ab1 == a & Ab2 == "Input" & Strain == s)
    Double = which(Ab1 == a & Ab2 == a & Strain == s)
    if( length(Input) == 1 && length(Double) == 1) 
    {
        NI = N.uniq[[Input]] / Summary["Input", "Mix.yield"]
        ND = N.uniq[[Double]] / Summary[a, "Mix.yield"]
        NI = NI/LoadingRatio[["Input"]]
        ND = ND/LoadingRatio[[a]]
        yield[[s]] = ND/NI
        r = ComputeTopRatio(Nucs[Input,], Nucs[Double,])
        r = r * Summary["Input", "Mix.yield"]/ Summary[a, "Mix.yield"]
        r = r * LoadingRatio[["Input"]] / LoadingRatio[[a]]
        ryield[[s]] = r
    }
  }
  print(a)
  print(unlist(yield))
  print(unlist(ryield))
  
  if( length(ryield) > 0)
    AbYield[[a]] = median(unlist(ryield),na.rm = TRUE)
}
print(AbYield)

SYield = list( WT = 1)
for( s in levels(Strain) )
  if( s != "WT" ) {
    ryield = list()
    for( a in levels(Ab1) )
      for( b in levels(Ab2))
      {
        WT = which(Ab1 == a & Ab2 == b & Strain == "WT")
        S = which(Ab1 == a & Ab2 == b & Strain == s)
        if( length(WT) == 1 && length(S) == 1) 
        {
          r = ComputeTopRatio(Nucs[WT,], Nucs[S,])
          ryield[[paste(a,b)]] = r
        }
      }
    print(s)
    print(unlist(ryield))
    if( length(ryield) > 0)
      SYield[[s]] = median(unlist(ryield), na.rm = TRUE)
  }

print(SYield)

TotalYield = array(dim=c(length(Names)))
dimnames(TotalYield)[[1]] = Names
for( i in 1:length(Names)) {
  seq = Summary[Ab2[[i]], "Mix.yield"]
  a1 = AbYield[[Ab1[[i]]]]
  a2 = AbYield[[Ab2[[i]]]]
  s = SYield[[Strain[[i]]]]
  TotalYield[[i]] = seq * s* a1*a2 * LoadingRatio[[Ab2[[i]]]]
}

print(TotalYield)

saveRDS(TotalYield,paste0(DataDir,"/TotalYield.rdata"))
