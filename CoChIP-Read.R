#
# Fix chr names 
#
FixChrName <- function( x ) {
  if( grepl("MT", x) ) x <- "M"; 
  
  if( !grepl("chr", x) ) x <- paste("chr",x,sep="")
  
  return(x)  
}

FixChrNames <- function( rname ) {
  ls = levels(rname)
  return( factor(as.vector(rname), ls, labels = as.character(sapply(ls,FixChrName)) ))
}
#
# Transform our data structure to genomic ranges data structure
#
library(GenomicRanges)

BuildGRanges <- function( b, Is  ) {
  Rs <- GRanges( seqnames = FixChrNames(b$rname[Is]), 
                 ranges=IRanges(start = b$mpos[Is], end = b$pos[Is]),
                 strand = b$strand[Is]
  )
  
}

##
ReadUniBAMFile <- function( BAMFile, param, extend = 150 ) {
  print( BAMFile )
  indexBam(BAMFile) # index reads in bam file
  
  ## Read in sorted indexed bam file
  # bam <- scanBam(BAMFile, param=param)
  
  ## Filter based on length
  #len = bam[[1]]$pos - bam[[1]]$mpos
  #  legalread = (len > 50) & (len < 1000) & !is.na(bam[[1]]$pos) & !is.na(bam[[1]]$mpos)
  #legalread = !is.na(bam[[1]]$pos) & !is.na(bam[[1]]$mpos)
  
  # hist(len, main=BAMFile, breaks="FD")
  GR <- granges(readGAlignments(BAMFile,param=param))
  l <- seqlevels(GR)
  fl <- as.character(sapply(l, FixChrName))
  sn <- factor( as.vector(seqnames(GR)), l, labels=fl)
  seqlevels(GR) <- fl
  seqnames(GR) <- sn
  GR <- resize(GR, width=extend,fix="start")
  print(paste0(BAMFile, ": ", length(GR), " total frag"))
  dat <- list( BAMFile = BAMFile, GR = GR )
  
}

ReadPairedBAMFile <- function( BAMFile, param ) {
  print( BAMFile )
  if( !file.exists(paste0(BAMFile,".bai")) )
    indexBam(BAMFile) # index reads in bam file
  
  ## Read in sorted indexed bam file
  # bam <- scanBam(BAMFile, param=param)
  
  ## Filter based on length
  #len = bam[[1]]$pos - bam[[1]]$mpos
  #  legalread = (len > 50) & (len < 1000) & !is.na(bam[[1]]$pos) & !is.na(bam[[1]]$mpos)
  #legalread = !is.na(bam[[1]]$pos) & !is.na(bam[[1]]$mpos)
  
  # hist(len, main=BAMFile, breaks="FD")
  GR <- granges(readGAlignmentPairs(BAMFile,param=param))
  l <- seqlevels(GR)
  fl <- as.character(sapply(l, FixChrName))
  sn <- factor( as.vector(seqnames(GR)), l, labels=fl)
  seqlevels(GR) <- fl
  seqnames(GR) <- sn
  print(paste0(BAMFile, ": ", length(GR), " total frag"))
  dat <- list( BAMFile = BAMFile, GR = GR )
}

getBAMParam <- function () {
  what <- c("pos", "strand", "rname","seq", "mpos")
  flag <- scanBamFlag(isUnmappedQuery = FALSE, 
                      isPaired = TRUE,
 #                     isSecondaryAlignment = FALSE,
                      isNotPassingQualityControls = FALSE)
  
  ScanBamParam(what = what, flag = flag) 
}


ReadBAMFiles <- function( fs, PairedEnd = TRUE ) {

  param <- getBAMParam()
  
  if( PairedEnd )
    rf <- function( fn ) { ReadPairedBAMFile( fn, param) }
  else
    rf <- function( fn ) { ReadUniBAMFile( fn, param) }
  
  rs <- lapply(fs,rf)
  
  reads = rs
}

ReadBAMDirectory <- function( BAMDir, PairedEnd = TRUE ) {
  fs <- list.files( BAMDir, "*.bam$", full.names = TRUE )
  ReadBAMFiles(fs, PairedEnd )
}
