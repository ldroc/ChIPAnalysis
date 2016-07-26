setwd("~/Dropbox/CoChIP/coChIP-Analysis")
source("CoChip-Functions.R")
source("NucAtlas.R")

DataDir = "~/Data/CoChIP/coChIP-Analysis/150909"
BamDir = paste0(DataDir,"/MergeBAM")

## prepere SGD/nuc structures

GetSGDandNucs()

mod.nucs = ReadHistoneModeAtlas("Weiner-HisMod.csv", nucs, c("H3K4me3", "H3K36me3"))


ComputeExpectChIP <- function( mnucs, mod, sizes, probs) {
  m = mcols(mnucs)
  ## argh - different type casts to get this to work
  w = as.vector(exp(as.matrix(m["Input"]) + as.matrix(m[mod])))
          
  SizeCoverage <- function( s, p )
  {
    d = s - 150
    tmp <- resize(mnucs, width = 150+2*d, fix = "center" )
    coverage(tmp,weight = w) * p
  }
  l = mapply(SizeCoverage, sizes, probs)
#  c = sqrt(Reduce('+',lapply(l,function(x) x*x)))
  c = Reduce('+', l)
}

RleListtoGranges <- function(rr) {
  cc <- unlist(lapply(rr,nrun))
  cc <- Rle(names(cc), cc)
  ir <- IRanges(start = unlist(lapply(rr, start)),
                width = unlist(lapply(rr, width)) )
  vs <- unlist(lapply(rr,runValue))
  GRanges(seqnames = cc, ranges = ir, mcols = vs)
}

GRtoDataTrack <- function( gr, focus, title) {
  gr.focus = gr[queryHits(findOverlaps(gr,focus))]
  DataTrack(RleListtoGranges(coverage(gr.focus)), type = "histogram",name = title) 
}

sizes = c( 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400 )
probs = c( .8, 1, 1.5, 2, 1.5, 1.2, 1, .8, .6, 0.5, 0.4, 0.3) 

focus <- GRanges(c("chrIII"), IRanges(31250,52750) ) 
#focus <- GRanges(c("chrIII"), IRanges(31250,32000) ) 
  focus.nucs <- mod.nucs[queryHits(findOverlaps(mod.nucs, focus))]
  
  H3K4.gr = dat[[1]]$GR
  len = Rle(sort(width(H3K4.gr)))
  K4sizes = seq(140,500,5)
  K4probs = unlist(lapply(sizes, function(x) sum(runLength(len)[runValue(len) >= x & runValue(len) < x+5])))
  K4probs = probs/sum(probs)

  erle = ComputeExpectChIP(focus.nucs, "H3K4me3", K4sizes, K4probs)
  eGR = RleListtoGranges(erle)
  PK4track <- DataTrack(RleListtoGranges(erle),type="histogram",name="Predict H3K4me3")

  H3K36.gr = dat[[4]]$GR
  len = Rle(sort(width(H3K36.gr)))
  K36sizes = seq(140,500,5)
  K36probs = unlist(lapply(sizes, function(x) sum(runLength(len)[runValue(len) >= x & runValue(len) < x+5])))
  K36probs = probs/sum(probs)
  erle = ComputeExpectChIP(focus.nucs, "H3K36me3", K36sizes, K36probs)
  eGR = RleListtoGranges(erle)
  PK36track <- DataTrack(RleListtoGranges(erle),type="histogram",name="Predict H3K36me3")
  
  K4K36.gr = dat[[6]]$GR
  AK4K36track <- GRtoDataTrack(K4K36.gr, focus,"K4me3-K36me3") 
  AK4track <- GRtoDataTrack(H3K4.gr, focus, "K4me3-Input") 
  AK36track <- GRtoDataTrack(H3K36.gr, focus, "K36me3-Input") 
  
    
  genestrack <- AnnotationTrack(SGD$Genes)
  gtrack <- GenomeAxisTrack()
  tmp <- exp(focus.nucs$Input)
  Inputtrack <- DataTrack(range=focus.nucs,data=tmp,type="histogram",name="Input")
  tmp <- exp(focus.nucs$H3K4me3 + focus.nucs$Input)
  K4track <- DataTrack(range=focus.nucs,data=tmp,type="histogram",name="K4me3")
  tmp <- exp(focus.nucs$H3K36me3 + focus.nucs$Input)
  K36track <- DataTrack(range=focus.nucs,data=tmp,type="histogram",name="K36me3")
  plotTracks(list(Inputtrack,K4track, PK4track, AK4track,
                  K36track, PK36track, AK36track, 
                  AK4K36track,genestrack), which=focus,from=31250,to=52750)
  plot(coverage(focus.nucs,weight=exp(focus.nucs$Input+focus.nucs$H3K4me3))$chrIII,type="l")
