library(rtracklayer)

NETDataDir = "~/Dropbox/CoChIP/NET-seq"

if( !exists("NETWigs"))
  NETWigs = list()

ReadNETData <- function() {
for( s in c("WT", "set2"))
  for( dir in c("plus","minus") ){
    Name = paste0(s, "_", dir)
    fn = paste0(NETDataDir,"/", Name, ".wig")
    ww = import.wig(fn)
    genome(ww) = genome(SGD$Genes)
    seqinfo(ww) = seqinfo(SGD$Genes)
    NETWigs[[Name]] <<- ww
  }
}


NetSeqSignal <- function(wig, Regions) {
  olaps <- findOverlaps(wig, Regions)
  df <- DataFrame(Nuc=subjectHits(olaps),value.value=score(wig)[queryHits(olaps)])
  scores <- aggregate(df, by=list(df$Nuc),sum)
  Scores = vector(length = length(Regions))
  Scores[] = 0
  #  names(Scores) = 1:length(NucRegions)
  Scores[scores$Group.1] = scores$value.value
  Scores
}

StrandedNetSeqSignal <- function( wigPos, wigNeg, Regions ){
  strands = as.vector(strand(Regions))
  ScorePos = NetSeqSignal(wigPos, Regions[strands == "+"])
  ScoreNeg = NetSeqSignal(wigNeg, Regions[strands == "-"])
  Score = rep(0,length(Regions))
  Score[strands == "+"] = ScorePos
  Score[strands == "-"] = ScoreNeg
  Score
}

GenerateNETTrans <- function() {
  GRs = list()
  for( n in names(NETWigs)) {
    ww = NETWigs[[n]]
    t = quantile(ww$score[1:20000], .5)
    print(paste(n,"t =",t))
    ww = ww[ww$score > t]
    cv = coverage(resize(ww,width = 75,fix = "center"))
    ranges = lapply(levels(seqnames(ww)), function(chr) { 
      sl = slice(cv[[chr]], lower=1); 
      ranges(sl[width(sl) > 500]) 
    })
    names(ranges) = levels(seqnames(ww))
    seqs = unlist(sapply(names(ranges), function(n) rep(n, length(ranges[[n]]))))
    GRs[[n]] <- GRanges(seqinfo = seqinfo(SGD$Genes),
                        ranges = IRanges(unlist(sapply(ranges,start)), 
                                         unlist(sapply(ranges,end))),
                        seqnames = seqs
    )
  }
  
  Set2Plus = GRs[["set2_plus"]][countOverlaps(GRs[["set2_plus"]],GRs[["WT_plus"]])==0]
  Set2Minus =  GRs[["set2_minus"]][countOverlaps(GRs[["set2_minus"]],GRs[["WT_minus"]])==0]
  strand(Set2Plus) = "+"
  strand(Set2Minus) = "-"
  Set2Transcripts = sort(c(Set2Plus, Set2Minus),ignore.strand=T)
  
  saveRDS(Set2Transcripts,paste0(NETDataDir,"/set2Trans.rdata"))
  export(Set2Transcripts,"~/Dropbox/CoChIP/160124/RPD3/Tracks/set2-transcripts.bed")
}