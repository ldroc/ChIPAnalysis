#library(EBImage)
library(Rsamtools);
library(ggplot2)
#library(gplots)
library(RColorBrewer)
library(corrplot)
library(MASS)
library(tools)
library(GenomicAlignments)
library(rtracklayer)
#library(Gviz)
#library(ggbio)

source("MixPoisson.R")
source("DensityScatter.R")
source("CoChIP-QC.R")
source("CoChIP-Read.R")
source("MetaGene.R")
source("SizeCoChIP.R")
source("V-Plot.R")
#source("coChIP-GLM.R")
#source("CoChIP-simulator.R")

ExprOrder <- function( gr ) {
  ord = order( gr$expr , na.last = TRUE)
  ord[1:(length(ord)[1] - length(which(is.na(gr$expr))))]
}

LengthOrder <- function( gr, genes ) {
  ord = order( width(genes) )
  neword = match(genes[ord]$acc,gr$acc)
  neword[ !is.na(neword) ]
}

ccParams <- function() {
  list( 
    Save = TRUE,
    DataDir = NULL,
    nuc = FALSE,
    cov = FALSE,
    cen = FALSE,
    meta = FALSE,
    doVPlot = FALSE,
    reuseSavedData = TRUE,
    MaxFragLen = 600,
    MinFragLen = 50,
    BAMParam = getBAMParam(),
    PairedEnd = TRUE,
    CenWidth = 50,
    Verbose = TRUE,
    TSSRegions = NULL,
    TTSRegions = NULL,
    GeneRegions = NULL,
    PlusOneRegions = NULL,
    NucRegions = NULL,
    NucRegionWidth = 50
  )
}

ccBuildFN <- function(Name, param, suff = ".rdata" ) {
  if( is.null(param$DataDir) )
    f = ""
  else
    f = paste0(param$DataDir,"/")
  paste0(f, Name, suff)
}

ExactNucCoverage <- function(NucRegions,GR) {
  NarrowNuc = resize(NucRegions,width=50,fix="center")
  olp = findOverlaps(NarrowNuc, GR, type="within")
  counts = as.table(t(olp))
  olp = findOverlaps(NarrowNuc, GR, type="any")
  counts.any = as.table(t(olp))
  MonoNucs = GR[counts == 1 & counts.any == 1]
  Mono = countOverlaps(NucRegions,MonoNucs)
  DiNucs = GR[counts == 2 & counts.any == 2]
  DiNucIndex = findOverlaps( DiNucs,NucRegions, select = "first")
  df = as.data.frame(table(DiNucIndex))
  Di = vector(length=length(NucRegions))
  Di[] = 0
  # there should be a better way to do this...
  Di[as.integer(as.character(df$DiNucIndex))]=df$Freq
  list( Mono = Mono, Di = Di)
}

ccProcessFile <- function( filename = NULL,
                           dat = NULL,
                           param = ccParams(),
                           Force = FALSE,
                           Change = FALSE
                    ) 
{
  Change = FALSE
  if( is.null( filename ) && is.null( dat ))
  {
    print("Need one of filename or dat be assigned!")
    return(NULL)
  }
  if( is.null( dat ) )
  {
    Name = BaseFileName(filename)
  } else {
    Name = dat$Name
    if( param$Verbose ) print(paste(Name, ": Processing precomputed data"))
  }
  fn = ccBuildFN(Name, param )
  if( is.null(dat) ) {
    if( param$reuseSavedData && file.exists(fn)) {
      if( param$Verbose ) print(paste(Name, ": Reading precomputed data"))
      dat <- readRDS(fn)
    } else {
        dat = list(Name = Name, 
                 GR = NULL, 
                 UniqGR = NULL,
                 nuc = NULL, 
                 cov = NULL, 
                 cen = NULL,
                 meta = NULL, 
                 vplot = NULL )
      Change = TRUE
    }
  }
    
  # Get BAM reads into a GenomicRanges object
  if( is.null( dat$GR) ) {
    if( param$Verbose ) print(paste(Name, ": Reading BAM file", filename))
    if( !file.exists(filename) )
      filename = paste0(filename, ".bam")
    if( param$PairedEnd )
      rs = ReadPairedBAMFile( filename, param$BAMparam)
    else
      rs <- ReadUniBAMFile( filename, param$BAMparam)
    dat$GR = rs$GR
    strand(dat$GR) = "*"
    Change = TRUE
  }

  if( Force ) {
    dat$UniqGR = NULL
    dat$cov = NULL
    dat$nuc = NULL
    dat$cen = NULL
    dat$meta = NULL
  }
  
  if( is.null(dat$UniqGR) ) {
    if( param$Verbose ) print(paste(Name, ": Building unique read set"))
    strand(dat$GR) = "*"
    I = width(dat$GR) < param$MaxFragLen & width(dat$GR) > param$MinFragLen
    dat$UniqGR = unique(dat$GR[I])
    dat$cov = NULL
    dat$nuc = NULL
    dat$cen = NULL
    Change = TRUE
  }
  
  # Raw coverage
  if( param$cov && is.null(dat$cov) ) {
    if( param$Verbose ) print(paste(Name, ": Computing coverage"))
    dat$cov = coverage(dat$UniqGR)
    Change = TRUE
  }
  
  # Center coverage
  if( param$cen && is.null(dat$cen) ) {
    if( param$Verbose ) print(paste(Name, ": Computing center coverage"))
    if( param$CenWidth > 0 ) {
      dat$cen = coverage(resize(dat$UniqGR,width=param$CenWidth,fix="center"))
    } else
      dat$cen = coverage(dat$UniqGR)
    
    Change = TRUE
  }
  
  if( param$nuc && is.null(dat$nuc) ) {
    print(paste(Name,"computing nucleosome coverage"))
    dat$nuc = unlist(countOverlaps(resize(param$NucRegions, width=param$NucRegionWidth, fix="center"), dat$UniqGR))
    Change = TRUE
  }
  
  if( param$meta && is.null(dat$meta) ) {
    if( is.null(dat$cen) )
      print(paste("Missing center coverage for file",Name))
    else {
      if( param$Verbose ) print(paste(Name, ": Computing gene centered coverage"))
      TSS <- CollectRegions(dat$cen, param$TSSRegions)
      TTS <- CollectRegions(dat$cen, param$TTSRegions)
      Genes <- CollectRegions(dat$cen, param$GeneRegions, width = 5000)
      PlusOne <- CollectRegions(dat$cen, param$PlusOneRegions)
      dat$meta = list( TSS= TSS, TTS=TTS, Genes=Genes, PlusOne = PlusOne )
      Change = TRUE
    }
  }
    
  # we are done!
  if( Change && param$Save ) {
    if( param$Verbose ) print(paste(Name, ": Saving data"))
    saveRDS(dat,fn)
  }
  
  return( dat )
}

ccDoQC <- function( dat, QCDir = NULL ) {
  if( is.null(QCDir) )
    fname = dat$Name
  else
    fname = paste0(QCDir,"/",dat$Name)
  fname = paste0(fname, ".png")
  QCGRanges( dat$Name, dat$GR, fname )
}

ccDoMeta <- function( dat, params, MetaDir = NULL, Qs = SGD$ExprQuan[5:1], UseGenes = FALSE, NameAdd = "" ) {
  # make sure we have the data we need
  if( is.null( dat$meta) ) {
    p = ccParams()
    p$meta = TRUE
    p$Save = FALSE # if you wanted to save this computation, should have done it beforehand
    dat = ccProcessFile( dat = dat, param = p )
  }
  
  TSS <- AvgRegionsSubgroups(dat$meta[[1]], params$TSSRegions, Qs) 
  TTS <- AvgRegionsSubgroups(dat$meta[[2]], params$TTSRegions, Qs) 
  Genes <- AvgRegionsSubgroups(dat$meta[[3]], params$GeneRegions,  Qs) 
  Plus1 <- AvgRegionsSubgroups(dat$meta[[4]], params$PlusOneRegions,  Qs) 

  if( UseGenes ) {
    PlotMultiCoverage( list(dat$Name), list(Genes), list(TTS), path = MetaDir, PDF = T, NameAdd = NameAdd )
  } else
    PlotMultiCoverage( list(dat$Name), list(Plus1), list(TTS), path = MetaDir, PDF = T , NameAdd = NameAdd)
}

ccPlotAlignedMatrix <- function( Dat, params, order = NULL, group = NULL, type = "TSS" ) {
  if( type == "TSS") {
    GR = params$TSSRegions
    Meta = Dat$meta[[1]]
    offset = 500
  } else
    if( type == "TTS" || type == "PostTTS") {
      GR = params$TTSRegions
      Meta = Dat$meta[[2]]
      offset = 1000
    } else
      if( type == "Plus1") {
        GR = params$PlusOneRegions
        Meta = Dat$meta[[4]]
        offset = 500
      } else
        if( type == "Genes") {
          GR = params$GeneRegions
          offset = 500
          Meta = Dat$meta[[3]]
        } else
          abort()
  
  Index = 1:length(GR)
  if( !is.null(order) )
    Index = order
  
  if( !is.null(group)) 
    Index = Index[is.element(GR$acc[Index], group)]
  
    
  print(length(Index))
  p = PlotHeatMap(t(Meta[Index,]), base = type, offset = offset, ylab = "") #, title = Dat$Name)
  return(p)
}

ccCombineDat <- function( dat1, dat2, Name = NULL ) {
  if( is.null(Name) )
    Name = dat1$Name
  
  dat = list(Name = Name, 
             GR = NULL, 
             nuc = NULL, 
             cov = NULL, 
             cen = NULL,
             meta = NULL, 
             vplot = NULL )
  if(!is.null(dat1$GR)) {
    dat$GR = c(dat1$GR, dat2$GR)
  }
  if(!is.null(dat1$nuc) && !is.null(dat2$nuc)) {
    dat$nuc = dat1$nuc + dat2$nuc
  }
  if(!is.null(dat1$cov) && !is.null(dat2$cov)) {
    dat$cov = dat1$cov + dat2$cov
  }
  if(!is.null(dat1$cen) && !is.null(dat2$cen)) {
    dat$cen = dat1$cen + dat2$cen
  }
  if(!is.null(dat1$meta) && !is.null(dat2$meta)) {
    dat$meta = mapply("+", dat1$meta, dat2$meta)
  }
  return(dat)
}

ccExportGeneSheet <- function( dat, param = ccParams(), Dir = NULL, type = "Genes", wSize = 10 )
{
  if( is.null(Dir) )
    fname = dat$Name
  else
    fname = paste0(Dir,"/",dat$Name)
  fname = paste0(fname, ".csv")
  
  if( is.null( dat$meta) ) {
    p = ccParams()
    p$meta = TRUE
    p$Save = FALSE # if you wanted to save this computation, should have done it beforehand
    dat = ccProcessFile( dat = dat, param = p )
  }
  if( type == "TSS") {
    GR = params$TSSRegions
    Meta = dat$meta[[1]]
    offset = 500
  } else
    if( type == "TTS" || type == "PostTTS") {
      GR = params$TTSRegions
      Meta = dat$meta[[2]]
      offset = 1000
    } else
      if( type == "Plus1") {
        GR = params$PlusOneRegions
        Meta = dat$meta[[4]]
        offset = 500
      } else
        if( type == "Genes") {
          GR = params$GeneRegions
          offset = 500
          Meta = dat$meta[[3]]
        } else
          abort()
  
  
  N = dim(Meta)[[1]]
  M = dim(Meta)[[2]]
  L = seq( 1, M, by = wSize)
  A = do.call(cbind, lapply(1:(length(L)-1), function(i) rowMeans(Meta[,L[[i]]:(L[[i+1]]-1)])))
  rownames(A) = GR$acc
  colnames(A) = as.character(L[1:(length(L)-1)] - offset -1)
  df = cbind(data.frame( YName = GR$acc, Name = GR$name, Description = GR$desc), as.data.frame(A, optional = FALSE))
  write.csv(x = df, file = fname, row.names = FALSE, quote = FALSE)
}

ccExportTrack <- function( dat, params, Tiles, Dir = params$DataDir){
  print("Exporting Track")
  temp = Tiles
  score(temp) = countOverlaps(Tiles,resize(dat$UniqGR,width=params$CenWidth,fix="center"))
  params.DataDir = Dir
  fname = ccBuildFN(dat$Name,params, suff = ".bw")
  print(fname)
  export( temp, fname, format="BigWig" )  
}
