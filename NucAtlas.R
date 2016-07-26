## 
## Read nucleosome atlas

GetSGDandNucs <- function() {
  if(file.exists("atlas.rdata")) {
    print("Loading Atlas")
    load(file="atlas.rdata",.GlobalEnv)
  } else {
    print("Read Atlas")
    ## SGD + TSS/TTS annotations
    SGD <<- ReadGeneAtlas("SGD_TSS.tab")
    ## read nucs
    NucRegions <<- ReadNucAtlas("Weiner-NucAtlas.csv", SGD )
    Nucs.cov <<- coverage(NucRegions)
    Plus1Nucs <<- GetNucAlignedRegions(NucRegions, SGD)
    GeneDirections <<- GetGeneDirection("gene-dir-annotation.csv")
    save(SGD, NucRegions, Nucs.cov, Plus1Nucs, GeneDirections, file="atlas.rdata")
  } 
}

ConvertFactor <- function( f, l ) {
  
}

ReadNucAtlas <- function( fn, SGD ) {
  df <- read.csv(fn, header = TRUE, stringsAsFactors = FALSE)

  lchrnames = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", 
               "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", 
               "chrXIV", "chrXV", "chrXVI" )
  # converge to GRanges
  
  chr = factor(df$chr, labels = lchrnames)
  acc = factor(df$acc, levels =levels(SGD$Genes$acc))
  rs <- GRanges( seqnames = chr,
                 ranges = IRanges( df$center - 75, df$center + 75),
                 acc = acc,
                 gene_pos = df$gene_pos,
                 seqinfo = Seqinfo(genome="sacCer3")
                 )
  
}



GetNucAlignedRegions <- function( nucs, SGD ) {
  gr <- nucs[nucs$gene_pos == "1"]
  getStrand <- function( a ) {
    I = (SGD$Genes$acc == a)
    I = which(I)
    if( length(I) > 0 )
      return( as.factor(strand(SGD$Genes)[I[1]]))
    else
      return(NA)
  }
  strand(gr) <- unlist(lapply(gr$acc, getStrand ))
  gr <- resize(gr, 2000, fix="center" )
  gr <- resize(gr, 1500, fix="end" )
}

ReadGeneAtlas <- function( fn ) {
  df <- read.csv(fn, sep="\t")
  
  chrnames = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", 
               "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", 
               "chrXIV", "chrXV", "chrXVI", "chrM" )
  lev = levels(df$chr)
  df$chr <- factor(df$chr, levels = lev, labels = chrnames)
  strand <- factor(replicate(length(df$ACC), "+"), levels=c("+", "-", "*"))
  strand[df$start > df$end] = "-"
  ustart = mapply(min, df$start, df$end)
  uend = mapply(max, df$start, df$end)
  seqinfo = Seqinfo(genome="sacCer3")
  I = df$Name == ""
  Names = as.character(df$Name)
  Names[I]= as.character(df$ACC[I])
  Desc = as.character(df$Description)
  gs <- GRanges( seqnames = df$chr,
                 ranges = IRanges(ustart, uend),
                 strand = strand,
                 acc = df$ACC,
                 name = Names,
                 expr = df$Expression,
                 desc = Desc,
                 seqinfo = seqinfo )
  
  GeneRegion = resize( resize( gs, width=4500, fix="start"),
                      width = 5000, fix="end" )
  
  Legal = !is.nan(df$TSS) & !is.nan(df$TTS)
  tss = df$TSS
#  tss[is.nan(tss)] = df$start[is.nan(tss)]
  tts = df$TTS
#  tts[is.nan(tts)] = df$end[is.nan(tts)]
  utss = mapply(min, tss, tts)
  utts = mapply(max, tss, tts)
  
  ## Fix situations where TSS/TTS are missing or wrong, 
  ## go with gene annotations in these cases
  B1 = utts < uend
  B2 = utss > ustart
  B = (B1 | B2) 
  B[is.na(B)] = TRUE
#  utss[B] = ustart[B]
#  utts[B] = uend[ B ]
  utss[B] = NA
  utts[B] = NA
  I = !B
  
  ts <- GRanges( seqnames = df$chr[I],
                 ranges = IRanges(utss[I],utts[I]),
                 strand = strand[I],
                 acc = df$ACC[I],
                 name = Names[I],
                 desc = Desc[I],
                 expr = df$Expression[I],
                 seqinfo = seqinfo )
    
  TSSRegion = resize(resize(ts,width=1000, fix="start"),
                       width=1500, fix="end" )
  TTSRegion = resize(resize(ts,width=1000, fix="end"),
                       width=1500, fix="start" )
  Plus = strand(ts) == "+"
  Neg = strand(ts) == "-"
  qs = quantile(ts$expr, probs = c(0,0.2,0.4,0.6,0.8,1), na.rm=TRUE)
  print(qs)
  ExprQuan = lapply( 1:5, 
                     function (i) ts$acc[ts$expr >= qs[i] & 
                                         ts$expr < qs[i+1] &
                                         !is.nan(ts$expr) ] )
  names(ExprQuan) = names(qs)[2:length(qs)]
  print(lapply(ExprQuan,length))
  return(list( Genes = gs, Transcripts = ts, 
               TSSRegion = TSSRegion, TTSRegion = TTSRegion, 
               GeneRegion = GeneRegion,
               Plus = Plus, Neg = Neg, 
               ExprQuan = ExprQuan))
}

AccToGeneName <- function( L, Genome ) {
  acc = Genome$Genes$acc
  name = Genome$Genes$name
  I = (name == "")
  NewN = as.character(name)
  NewN[I] = as.character(acc[I])
  Acc = as.character(acc)
  
  NewN[sapply(as.character(L), function(a) which(Acc == a)[[1]])]
}

ReadHistoneModeAtlas <- function( fn, nucs, modlist) {
  df <- read.csv(fn, header = TRUE, stringsAsFactors = FALSE)
  Inputs <- colnames(df)[grep("Input", colnames(df))]
  Input <- rowMeans(df[Inputs])
  Input <- Input - mean(Input)
  mcols(nucs, use.names = FALSE) <- cbind( mcols(nucs),  Input, df[modlist])
  nucs
}

GetGeneDirection <- function(fn) {
  df <- read.csv(fn, header = TRUE, stringsAsFactors = FALSE)
  list (
    Promoter.Continue = df$ORF[df$Promoter == 1],
    Promoter.Diverge = df$ORF[df$Promoter == -1],
    Promoter.Other = df$ORF[df$Promoter == 0],
    Terminator.Continue = df$ORF[df$Terminator == 1],
    Terminator.Diverge = df$ORF[df$Terminator == -1],
    Terminator.Other = df$ORF[df$Terminator == 0] )
  
}