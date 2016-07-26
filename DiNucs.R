WorkDir = "~/Dropbox/CoChIP/coChIP-Analysis"
setwd(WorkDir)


#DataDir =  "~/Data/CoChIP/HisMut-K4/HisMut"  
#DataDir =  "~/Data/CoChIP/160124/RPD3"  
DataDir = "~/Google Drive/CoChIPAnalysis"
FigureDir = "~/Dropbox/CoChIP/Analysis/TempFigures"
FigureDir = "~/Google Drive/CoChIPAnalysis/Figures"
TrackDir = "~/Google Drive/CoChIPAnalysis/DiNucTracks"

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

Sample = sample(1:length(DiNucID2),5000)
Sample = 1:length(DiNucID2)
s = "WT"
a1 = "K18ac"
a2 = "Input"

if(0) {
  for( a1 in TAb )
    for( a2 in TAb ) {
      title = paste0(s,"-",a1,"-",a2)
      p1 = DensityScatter(pmax(mono.Nucs[title, DiNucID1[Sample]],mono.Nucs[title, DiNucID2[Sample]]), 
                          di.Nucs[title,Sample], title = paste(title, "vs max", title), jitter = TRUE,
                          xlabel = "Mono", ylabel="DiNuc")
      p2 = DensityScatter(mono.Nucs[title, DiNucID1[Sample]]+mono.Nucs[title, DiNucID2[Sample]], 
                          di.Nucs[title,Sample],title = paste(title, "vs sum", title), jitter = TRUE,
                          xlabel = "Mono", ylabel="DiNuc")
      title1 = paste0(s,"-",a1,"-Input")
      title2 = paste0(s,"-",a2,"-Input")
      p3 = DensityScatter(pmax(mono.Nucs[title1, DiNucID1[Sample]],mono.Nucs[title1, DiNucID2[Sample]], 
                               mono.Nucs[title2, DiNucID1[Sample]],mono.Nucs[title2, DiNucID2[Sample]]),
                          di.Nucs[title,Sample],title = paste(title, "vs max individuals", title), jitter = TRUE,
                          xlabel = paste("max(", title1, ",", title2,")"), ylabel="DiNuc")
      p4 = DensityScatter(mono.Nucs[title1, DiNucID1[Sample]]+mono.Nucs[title1, DiNucID2[Sample]]+
                            mono.Nucs[title2, DiNucID1[Sample]]+mono.Nucs[title2, DiNucID2[Sample]],
                          di.Nucs[title,Sample],title = paste(title, "vs sum individuals", title), jitter = TRUE,
                          xlabel = paste("sum(", title1, ",", title2,")"), ylabel="DiNuc")
      fn = paste(FigureDir,"/di-vs-mono-",title,".png")
      png(fn, width=400, height=1200)
      multiplot(p1,p2,p3,p4)
      dev.off()
    }
}

library(rtracklayer)
BAMDir = "~/Data/CoChIP/160124/RPD3/BAM/"
Files = as.list(sapply(TAb, function(a) paste0("bar1-",a, "-",c(TAb, "Input"))))

if(0) {
  
  dats = lapply(Files, function(f) ccProcessFile(filename = paste0(BAMDir,f), param = params ))
  NarrowNuc = resize(NucRegions,width=50,fix="center")
  export(NucRegions, paste0(TrackDir,"/nuc-atlas.bed"))
  export(NarrowNuc, paste0(TrackDir,"/nuc-center-atlas.bed"))
  for( d in dats) {
    GR = d$UniqGR
    olp = findOverlaps(NarrowNuc, GR, type="within")
    counts = as.table(t(olp))
    olp = findOverlaps(NarrowNuc, GR, type="any")
    counts.any = as.table(t(olp))
    MonoNucs = GR[counts == 1 & counts.any == 1]
    DiNucs = GR[counts == 2 & counts.any == 2]
    NoNucs = GR[counts != counts.any | counts == 0 | counts > 2]
    export(MonoNucs, paste0(TrackDir,"/", d$Name, "-mono.bam"))
    export(DiNucs, paste0(TrackDir,"/", d$Name, "-di.bam"))
    export(NoNucs, paste0(TrackDir,"/", d$Name, "-none.bam"))
  }

  nucs.list = lapply(dats,function(d) ExactNucCoverage(NucRegions.select, d$UniqGR))
  Names = lapply(dats,function(d) d$Name)
  Names = sub("bar1-", "WT-", Names)
  M.Mono <- do.call(rbind,lapply(nucs.list, function(n) n$Mono))
  M.Di <- do.call(rbind,lapply(nucs.list, function(n) n$Di))
  rownames(M.Mono) = Names
  rownames(M.Di) = Names
}

# QC 
if(0) {
  p = params
  p$Save = TRUE
  p$meta = TRUE  
  p$cen = TRUE
  params.short = p
  params.long = p
  params.short$MaxFragLen = 180
  params.long$MinFragLen = 180
  
  if( !exists("dats")) {
    dats = list()
    for( f in Files )
      dats[[f]] = readRDS(paste0(DataDir,"/Data/",f,".rdata"))
  }
  
  for( d in dats ) {
    d.short = d
    d.short$Name = paste0(d$Name,"-short")
    d.short = ccProcessFile(dat=d.short,param=params.short, Force = TRUE)
    ccDoMeta(d.short, params, MetaDir = FigureDir, UseGenes = TRUE)
    ccDoMeta(d.short, params, MetaDir = FigureDir, Qs = LongGenes.Qs[5:1], UseGenes = TRUE, NameAdd = "-LongGenes" )
    d.long = d
    d.long$Name = paste0(d$Name,"-long")
    d.long = ccProcessFile(dat=d.long,param=params.long, Force = TRUE)
    ccDoMeta(d.long, params, MetaDir = FigureDir, UseGenes = TRUE)
    ccDoMeta(d.long, params, MetaDir = FigureDir, Qs = LongGenes.Qs[5:1], UseGenes = TRUE, NameAdd = "-LongGenes" )
  }
  for( d in dats ) 
    ccDoQC(d, FigureDir)
}

#
# tiled occupancy
#
if(0) {
  TileWidth = 20
  Tiles = tileGenome(seqinfo(SGD$Genes),tilewidth = TileWidth,cut.last.tile.in.chrom = TRUE)
  Shorts = list()
  Longs = list()
  Ratios = list()
  
  for( f in Files) {
      d = readRDS(paste0(DataDir,"/Data/",f,".rdata"))
      W = width(d$UniqGR)
      cov.short = unlist(countOverlaps(Tiles, d$UniqGR[W < 180]))
      cov.long = unlist(countOverlaps(Tiles, d$UniqGR[W > 180 & W < 400  ]))
      Shorts[[f]]  = cov.short
      Longs[[f]] = cov.long
      temp = Tiles
      score(temp) = cov.short
      export(temp, paste0(TrackDir,"/",f,"-short.bw"))
      score(temp) = cov.long
      export(temp, paste0(TrackDir,"/",f,"-long.bw"))
      Ratio = pmax(-20,pmin(log2(cov.long+.001) - log2(cov.short +.001),20))
      Ratio[is.na(Ratio)] = 0
      Ratios[[f]] = Ratio
      score(temp) = Ratio
      export(temp, paste0(TrackDir,"/",f,"-ratio.bw"))
  }
  
  CollectMeta <- function(W) {
    cov = coverage(Tiles,weight = W)
    meta.Gene = CollectRegions(cov, SGD$GeneRegion, width = 5000)
    meta.TTS = CollectRegions(cov, SGD$TTSRegion, width = 1500)
    list( Gene = meta.Gene, TTS = meta.TTS )
  }
  
  PlotMeta <- function(cov, fname) {
    meta = CollectMeta(cov)
    TTS.all <- AvgRegionsSubgroups(meta[["TTS"]], SGD$TTSRegion, SGD$ExprQuan[5:1]) 
    TTS.long <- AvgRegionsSubgroups(meta[["TTS"]], SGD$TTSRegion, LongGenes.Qs[5:1]) 
    Gene.all <- AvgRegionsSubgroups(meta[["Gene"]], SGD$GeneRegion, SGD$ExprQuan[5:1]) 
    Gene.long <- AvgRegionsSubgroups(meta[["Gene"]], SGD$GeneRegion, LongGenes.Qs[5:1]) 
    
    PlotMultiCoverage( list(fname, paste0(fname,"-longGenes")), 
                       list(Gene.all, Gene.long), 
                       list(TTS.all, TTS.long), 
                       path = FigureDir, PDF = T )
  }
  
  for( f in Files ) {
    PlotMeta(Shorts[[f]], paste0("tiles-",f,"-short"))
    PlotMeta(Longs[[f]], paste0("tiles-",f,"-long"))
    PlotMeta(Ratios[[f]], paste0("tiles-",f,"-ratio"))
  }
    
}

# length histogram per position
if( 0 ) {
  PlotHistogram <- function(d, Regions, fname, binwidth = 5) {
    W = width(d$UniqGR)
    df = data.frame( x = W, color=rep(0,length(W)))
    for( i in 1:length(Regions) ) {
      n = countOverlaps(d$UniqGR,Regions[[i]])
      df[n > 0, "color"] = as.character(i)
    }
    p = ggplot(df)
    p = p + scale_color_manual(values = rainbow(length(Regions)),
                           labels=names(Regions),
                             guide = "legend")
    p = p + geom_histogram(data=df,
                           binwidth = binwidth,
                           aes(x = x, y = ..density..), 
                           color="darkgray")
    
    for( i in 1:length(Regions) ) {
      p = p + geom_freqpoly(data=df[df$color == i,],
                            binwidth = binwidth,
                            aes(x = x, y = ..density.., color=color), 
                            show.legend = TRUE,size = .5)
    }
    p = p + labs(title=fname, x="width", y="frequency")
    p = p+theme(legend.position="right")
    ggsave(plot = p, filename = paste0(FigureDir,"/", fname, ".pdf"))
  }

  ComputeNucRatio  <- function( d, Regions ) {
    GR.short = d$UniqGR[ width(d$UniqGR) <= 200]
    rs = sapply(Regions, function(r) sum(countOverlaps(GR.short, r)>0)/sum(countOverlaps(d$UniqGR,r)>0) )
    rs = c(Background = sum(countOverlaps(GR.short, NucRegions.select)>0)/sum(countOverlaps(d$UniqGR,NucRegions.select)>0), rs)
  }
  NucPos = c(-1,1,2,3,4,5,6)
  Regions = lapply(NucPos, 
                   function(n) {
                     I = !is.na(NucRegions.select$gene_pos) & NucRegions.select$gene_pos == n
                     NucRegions.select[I]
                   })
  Regions = c(Regions,
              LongGeneBody = NucRegions.select[!is.na(NucRegions.select$gene_pos) & NucRegions.select$gene_pos > 6],
              NonGene = NucRegions.select[is.na(NucRegions.select$gene_pos)] )
  names(Regions) = c(paste("Nuc",NucPos),"LongGeneBody", "NoneGene")
  Nuc1Ratio = list()
  for( f in Files) {
    d = readRDS(paste0(DataDir,"/Data/",f,".rdata"))
    if( length(d$UniqGR) > 25000 ) {
      PlotHistogram(d,Regions, paste0("hist-",f)) 
      Nuc1Ratio[[f]] = ComputeNucRatio(d,Regions)
    }
  }
  write.csv(t(as.matrix(as.data.frame(Nuc1Ratio))), file = paste0(DataDir,"/Nuc1Ratios.csv"))
}