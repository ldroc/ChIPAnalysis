setwd("~/Dropbox/CoChIP/coChIP-Analysis")
source("CoChip-Functions.R")
source("NucAtlas.R")
source("Multiplot.R")
source("NucAttr.R")

DataDir = "~/Dropbox/CoChIP/Analysis"
FigureDir = "~/Dropbox/CoChIP/Analysis/TempFigures"

source("PairwiseModel.R")

source("NET-Seq.R")
NETDataDir = "~/Dropbox/CoChIP/NET-seq"
library(rtracklayer)

TranscriptTracks = list( xuts = "van_Dijk_2011_XUTs_V64.bed", 
                         cuts = "Neil_2009_class_I_CUTs_V64.bed")
#                         ncRNAs = "Neil_2009_class_I_ncRNAs_V64.bed" )

LoadCUTtracks <- function() {
  ProcessTrack <- function(t) {
    r = import(paste0(NETDataDir,"/",t))
    r$score = StrandedNetSeqSignal(NETWigs$WT_plus, NETWigs$WT_minus, r)
    r$acc = r$name
    r
  }
  Tracks <- lapply(TranscriptTracks, ProcessTrack)
  Transcripts <- SGD$Transcripts
  Transcripts$score = StrandedNetSeqSignal(NETWigs$WT_plus, NETWigs$WT_minus, Transcripts)
  Genes <- SGD$Genes[seqnames(SGD$Genes) != "chrM"]
  Genes$score = Genes$expr
  Tracks <- c(Tracks, list(Transcripts = Transcripts,  Genes = Genes) )
}

if(length(NETWigs) == 0)
  ReadNETData()

Tracks = LoadCUTtracks()

Thresholds = lapply(Tracks, function(t) quantile(t$score, c(1, .8, .6, .4, .2),na.rm = TRUE))

BuildSubGroups <- function( t, Ts) {
  Groups = lapply(1:(length(Ts)-1), function(i) t$acc[which(t$score <= Ts[[i]] & t$score > Ts[[i+1]])])
  names(Groups) = lapply(1:(length(Ts)-1), function(i) paste0(names(Ts)[[i+1]],"-",names(Ts)[[i]]))
  Groups
}

SubGroups = mapply(BuildSubGroups, Tracks, Thresholds,SIMPLIFY = FALSE)
TRegions = lapply(Tracks, function(t) resize(resize(t,width=2000, fix="start"),width=2500,fix="end"))
                  
NucCoverage <- function(N) {
  cov = coverage(NucRegions.select,weight=N)
}

BuildTrackMeta <- function( cov, title ) {
  meta = lapply(TRegions, function(t) CollectRegions(cov, t, width = 2500))
  As = mapply(AvgRegionsSubgroups, meta, Tracks, SubGroups, SIMPLIFY = FALSE)
  fn = paste0(FigureDir,"/cuts-meta-",title,".png")
  png(fn,height = 1024, width=600)
  parparam = par()
  n = length(Tracks)
  par(mfrow = c(n,1))
  for( i in 1:n) {
    main =   paste(names(Tracks)[[i]], title)
    PlotCovergeGroups(As[[i]], main=main, label="TSS", pos=500)
  }
  dev.off()
  suppressMessages(par( mfrow = parparam$mfrow))
}

BuildTrackNET <- function( ) {
  covPlus = coverage(NETWigs$WT_plus,weight = score(NETWigs$WT_plus))
  covNeg = coverage(NETWigs$WT_minus,weight = score(NETWigs$WT_minus))
  
  fixmeta <- function(t, mP, mN) {
    meta = mP
    strand = as.vector(strand(t))
    meta[strand == "-",] = mN[strand == "-",]
    meta
  }

  metaP = lapply(TRegions, function(t) CollectRegions(covPlus, t, width = 2500))
  metaN = lapply(TRegions, function(t) CollectRegions(covNeg, t, width = 2500))
  metaS = mapply(fixmeta, TRegions, metaP, metaN)
  metaAS = mapply(fixmeta, TRegions, metaN, metaP)
  
  AvgS = mapply(AvgRegionsSubgroups, metaS, Tracks, SubGroups, SIMPLIFY = FALSE)
  AvgAS = mapply(AvgRegionsSubgroups, metaAS, Tracks, SubGroups, SIMPLIFY = FALSE)
  As = list( Sense = AvgS, AntiSense = AvgAS)
  YMax = c( 5, 10, 15 )
  for( j in 1:length(As))
    for( y in YMax ) {
      title = paste0("NET-seq-",names(As)[[j]])
      fn = paste0(FigureDir,"/cuts-",title,"-",y,".png")
      print(fn)
      png(fn,height = 1024, width=600)
      parparam = par()
      n = length(Tracks)
      par(mfrow = c(n,1))
      for( i in 1:n) {
        main =   paste(names(Tracks)[[i]], title)
        PlotCovergeGroups(As[[j]][[i]], main=main, label="TSS", pos=500, ylim=c(0,y))
      } 
      dev.off()
      suppressMessages(par( mfrow = parparam$mfrow))
    }
}


s = "WT"
if(1) {
  for( i in 1:length(Tracks))
    BuildTrackMeta(coverage(Tracks[[i]]),names(Tracks)[[i]])
  
  Seen = list()
  for( s in c("WT", "2.WT", "2.eaf3", "2.set2"))
    for( a1 in TAb )
      for( a2 in TAb )
        if( a1 != a2 ) {
          print(paste(s,a1,a2))
          mod = ComputePairwiseExpectations(s,a1,a2)
          BuildTrackMeta(NucCoverage(mod$V.obs), paste0(s,"-", a1,"-", a2))
          if( !(a1 %in% names(Seen))) {
            BuildTrackMeta(NucCoverage(Nuc(s,a1)/mod$t1), paste0(s,"-", a1,"-", "Input"))
            Seen[[a1]] = 1
          }
        }
}

BuildTrackNET() 