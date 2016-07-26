library(Rsamtools);
library(ggplot2)
library(gplots)
library(RColorBrewer)
setwd("~/Dropbox/CoChIP/coChIP-Analysis")

x = source("NucAtlas.R")
source("CoChip-Functions.R")

## SGD + TSS/TTS annotations
SGD <- ReadGeneAtlas("SGD_TSS.tab")



GetRNASeqRegions <- function(ts, Len) {
    tsLen = resize(ts,width=Len, fix = "end")
    s = start(ranges(tsLen))
    e = end(ranges(tsLen))
    zz = mapply(max, start(ranges(ts)),start(ranges(tsLen)))
    I = as.vector(strand(ts) == '+')
    s[I] = zz[I]
    zz = mapply(min, end(ranges(ts)),end(ranges(tsLen)))
    I = as.vector(strand(ts) == '-')
    e[I] = zz[I]
    ranges(tsLen) = IRanges(s, e)
    return(tsLen)
}

Len = 750
tsGenes <- GetRNASeqRegions(SGD$Genes, Len)
tsTrans <- GetRNASeqRegions(SGD$Transcripts, Len)


rs <- ReadBAMDirectory("~/Data/RNA_seq/150730", PairedEnd = FALSE )

GetXprCount <- function( gr, ts ) {
    x <- countOverlaps(ts, gr)
    N <- length(gr) / 10**6
    x <- x / N
}

xs <- lapply( rs, function( r ) GetXprCount( r$GR, tsTrans) )
names(xs) <- lapply(rs, function( r) BaseFileName(r$BAMFile) )

xdat <- data.frame(  ACC = tsTrans$acc, GeneName = tsTrans$name, xs  )

write.csv(xdat, "~/Data/RNA_seq/150730.csv", row.names = FALSE)

eps = 1
RIplot <- function( x, y, xname, yname ) {
  I = (x+y)/2
  R = (x+eps)/(y+eps)
  DensityScatter(I, R, alpha=0.5, title=paste(xname, "vs", yname),
                 threshold=1, contour=FALSE, xlabel = "Intensity", ylabel="Ratio",
                 ylim=c(0.1,10),xlim=c(10,max(I)),
                 logscale = TRUE, horizontal=1, coordeq = FALSE, smooth=TRUE)
}


if( 0 ) {
frag = rs[[1]]$GR
x = countOverlaps(rs[[1]]$GR, tsTrans)
x1 = countOverlaps(rs[[1]]$GR, tsTransLong)
mfrag = frag[x > 1]
mfrag1 = frag[x1 > 1]
y = countOverlaps(tsTrans,rs[[1]]$GR, ignore.strand = FALSE)
y1 = countOverlaps(tsGenes,rs[[1]]$GR, ignore.strand = FALSE)
y3 = countOverlaps(tsTransLong,rs[[1]]$GR, ignore.strand = FALSE)

z = countOverlaps(tsTrans, mfrag, ignore.strand = FALSE)
z1 = countOverlaps(tsTransLong, mfrag1, ignore.strand = FALSE)

}