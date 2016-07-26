WorkDir = "~/Dropbox/CoChIP/coChIP-Analysis"
setwd(WorkDir)



DataDir = "~/Google Drive/CoChIPAnalysis"
FigureDir = "~/Google Drive/CoChIPAnalysis/K4Figures"
TrackDir = "~/Google Drive/CoChIPAnalysis/K4Tracks"

print("Initializing")
suppressMessages(source("CoChip-Functions.R"))
suppressMessages(source("NucAtlas.R"))
suppressMessages(source("DensityScatter.R"))
GetSGDandNucs()

source("PairwiseModel.R")

params = ccParams()
params$DataDir = paste0(DataDir,"/Data")
params$cov = FALSE
params$meta = FALSE
params$cen = FALSE
params$nuc = FALSE
params$TSSRegions = SGD$TSSRegion
params$TTSRegions = SGD$TTSRegion
params$GeneRegions = SGD$GeneRegion
params$PlusOneRegions = Plus1Nucs
params$NucRegions = NucRegions.select
params$MinFragLen = 0
params$MaxFragLen = 1000

##
## BuildDiNucAtlas
##
DiNucFN = paste0(WorkDir,"/DiNucRegions.rdata")


LongGenes = width(SGD$Genes) > 2000 & !is.na(width(SGD$Genes)) 
LongGenes.Expr = SGD$Genes[LongGenes]$expr
LongGenes.Expr[is.na(LongGenes.Expr)] = 0
qs = quantile(LongGenes.Expr, c(0,0.2,0.4,0.6,0.8,1), na.rm=TRUE)
qs[[6]] = Inf
LongGenes.TempQs = findInterval(LongGenes.Expr,qs)
#  LongGenes.Qs = lapply(qs[2:6], function(q) LongGenes)
LongGenes.Qs = list()
for( i in 1:5) {
  LongGenes.Qs[[i]] = SGD$Genes$acc[LongGenes][LongGenes.TempQs == i]
}
names(LongGenes.Qs) = c("20%", "40%", "60%", "80%", "100%")

if( !file.exists(DiNucFN)) {
  start = c()
  end = c()
  seq = c()
  names = c()
  for( chr in levels(seqnames(NucRegions.select)) )
  {
    cnames = which(seqnames(NucRegions.select) == chr)
    cNucs = NucRegions.select[seqnames(NucRegions.select) == chr]
    cstart = start(cNucs)
    cend = end(cNucs)
    N = length(cstart)
    if( N > 1 ) {
      cstart = cstart[2:N]
      cend=cend[1:(N-1)]
      cnames = cnames[1:(N-1)]
      I = (cend - cstart) < 350 
      cstart = cstart[I]
      cend = cend[I]
      cnames = cnames[I]
      N = length(cstart)
      start = c(start,cstart)
      end = c(end,cend)
      seq = c(seq,rep(chr,N))
      names = c(names,cnames)
    }
  }
  DiNucRegions = GRanges(seqnames = seq, 
                         ranges = IRanges(start = start, end = end, names=names),
                         seqinfo = seqinfo(SGD$Genes))
  
  saveRDS(DiNucRegions, DiNucFN)
} else
  DiNucRegions = readRDS(DiNucFN)

DiNucID1 = as.integer(names(DiNucRegions))
DiNucID2 = DiNucID1+1

DiNucDataFN = paste0(DataDir,"/DiNuc.rdata")
if( file.exists(DiNucDataFN) ) {
  ll= readRDS(DiNucDataFN)
  di.Nucs = ll$di
  mono.Nucs = ll$mono
} else {
  DataDirs = list( HisMut = "~/Data/CoChIP/HisMut-K4/HisMut/Data/",
                   RPD3 = "~/Data/CoChIP/160124/RPD3/DataLong/", 
                   K4Sym = "~/Data/CoChIP/160124/K4Sym/DataLong/" )

  ReadNucs <- function(f) {
    print(f)
    d = readRDS(paste0(dir,f))
    strand(d$GR) = "*"
    UniqGR = unique(d$GR)
    
    ExactNucCoverage(NucRegions.select, UniqGR)
  }
  
  i = 1
  di.Nucs.matrix = list()
  mono.Nucs.matrix = list()
  for( Expr in names(DataDirs) ) {
    dir = DataDirs[[Expr]]
    Files = list.files( dir, "*rdata$", full.names = F )
    Files =  c(Files[grep("WT", Files)],  Files[grep("bar1", Files)])
    nucs.list = lapply(Files,ReadNucs)
    M.Mono <- do.call(rbind,lapply(nucs.list, function(n) n$Mono))
    M.Di <- do.call(rbind,lapply(nucs.list, function(n) n$Di))
    rownames(M.Mono) = sub(".rdata","", Files)
    rownames(M.Mono) <- sub("bar1", "WT", rownames(M.Mono))
    rownames(M.Mono) <- gsub("H3K", "K", rownames(M.Mono))
    rownames(M.Mono) = paste0(i,".",rownames(M.Mono))
    rownames(M.Di) = rownames(M.Mono)
    i = i+1
    di.Nucs.matrix[[Expr]] = M.Di
    mono.Nucs.matrix[[Expr]] = M.Mono
  }
  di.Nucs = do.call(rbind,di.Nucs.matrix)
  di.Nucs = AggregateNucs(di.Nucs)
  di.Nucs = di.Nucs[,DiNucID1]
  colnames(di.Nucs) = colnames(Nucs)[DiNucID1]
  mono.Nucs = do.call(rbind, mono.Nucs.matrix)
  mono.Nucs = AggregateNucs(mono.Nucs)
  colnames(mono.Nucs) = colnames(Nucs)
  saveRDS(list(mono = mono.Nucs, di = di.Nucs), DiNucDataFN)
} 

TAb = c( "K4me1", "K4me2", "K4me3", "K4ac")
s = "3.WT"

orig.Nucs = Nucs
Nucs = mono.Nucs

if(0) {
  Nucs = mono.Nucs
  for( a1 in TAb )
    for( a2 in TAb) {
      print(paste("Check pair ",s,a1,a2))
      png(paste0(FigureDir, "/check-model-",s,"-",a1,"-",a2,".png"), width=400, height=1200)
      CheckPlots(s,a1,a2)
      dev.off()
    }
}

if(1) {
  Nucs = mono.Nucs
  for( a1 in TAb )
    for( a2 in c(TAb,"Input") ) {
        cov = coverage(NucRegions.select,weight = ModelDistribution(s, a1, a2))
        meta.Gene = CollectRegions(cov, SGD$GeneRegion, width = 5000)
        genes <- AvgRegionsSubgroups(meta.Gene, SGD$GeneRegion,  SGD$ExprQuan[5:1]) 
        pdf( paste0(FigureDir,"/meta-mono-",s,"-",a1,"-",a2,".pdf"))
        PlotCovergeGroups(lapply(genes, function(g) g[1:2500]))
        dev.off()
    }
  Nucs = mono.Nucs
  Nucs[,] = 0
  Nucs[,DiNucID1] = di.Nucs[,]
  Nucs[,DiNucID2] = Nucs[,DiNucID2] + di.Nucs[,]
  
  for( a1 in TAb )
    for( a2 in c(TAb,"Input") ) {
      cov = coverage(NucRegions.select,weight = ModelDistribution(s, a1, a2))
      meta.Gene = CollectRegions(cov, SGD$GeneRegion, width = 5000)
      genes <- AvgRegionsSubgroups(meta.Gene, SGD$GeneRegion,  SGD$ExprQuan[5:1]) 
      pdf( paste0(FigureDir,"/meta-di-",s,"-",a1,"-",a2,".pdf"))
      PlotCovergeGroups(lapply(genes, function(g) g[1:2500]))
      dev.off()
    }
}
