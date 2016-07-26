
## prepere SGD/nuc structures
if( !exists("SGD"))
  GetSGDandNucs()
source("Main-Functions.R")


Sample = 1:10000
Sample = 1:65032
Nuc <- function( s, a1, a2 = "Input") {
  Z = Nucs[paste0(s,"-",a1,"-",a2),] /Inputs
  #  Z = Nucs[paste0(s,"-",a1,"-",a2),]
  Z[Sample] #/(Seq.Yield[[a2]]*Ab.Spec[[a1]]*Ab.Spec[[a2]]*StrainRatio[[s]])
}
#TAb = c("K4me3", "K18ac", "K36me3", "K79me3", "K56ac")
TAb = c("K4me3", "K18ac", "K36me3", "K79me3")

AggregateNucs <- function( Nucs, 
                           TAb = c("K4me3", "K18ac", "K36me3", "K79me3", "K56ac"),
                           WT = "WT",
                           WTStrains = c("1.WT", "2.WT") ) {
  Names = rownames(Nucs)
  Strain <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[1]])))
  Ab1 <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[2]])))
  Ab2 <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[3]])))
  
  
  ## merge libraries
  ##
  MNucs <- matrix(0, nc = dim(Nucs)[[2]], nr = length(TAb) * (length(TAb)+1))
  colnames(MNucs) = colnames(Nucs)
  rownames(MNucs) = unlist(lapply( TAb, function(a1) paste0(WT,"-",a1,"-",c(TAb,"Input"))))
  for( a1 in TAb )
    for( a2 in c(TAb, "Input"))
    {
      r = paste0(a1,"-",a2)
      for( s in WTStrains )
        if( paste0(s,"-",r) %in% Names)
          MNucs[paste0(WT,"-",r),] = MNucs[paste0(WT,"-",r),] + Nucs[paste0(s,"-",r),]
    }
  
  Nucs = rbind( MNucs, Nucs)
}

if( !exists("Nucs")) {
  NucsFN = paste0(DataDir, "/Nucs.rdata")
  Nucs <- readRDS(NucsFN)
  
  if(0) {
  NucsLongFN = paste0(DataDir, "/Nucs-long.rdata")
  if( file.exists(NucsLongFN) ) {
    NucsLong <- readRDS(NucsLongFN)
  } else {
    K4SymNuc <- readRDS(paste0(DataDir,"/K4Sym-long-nuc.rdata"))
    HisMutNuc <- readRDS(paste0(DataDir,"/HisMut-long-nuc.rdata"))
    RPD3Nuc <- readRDS(paste0(DataDir,"/RPD3-long-nuc.rdata"))
    
    rownames(K4SymNuc) <- sub("bar1", "WT", rownames(K4SymNuc))
    rownames(RPD3Nuc) <- sub("bar1", "WT", rownames(RPD3Nuc))
    rownames(HisMutNuc) <- gsub("H3K", "K", rownames(HisMutNuc))
    
    rownames(HisMutNuc) <- paste0("1.", rownames(HisMutNuc))
    rownames(RPD3Nuc) <- paste0("2.", rownames(RPD3Nuc))
    rownames(K4SymNuc) <- paste0("3.", rownames(K4SymNuc))
    
    NucsLong <- rbind( HisMutNuc, RPD3Nuc, K4SymNuc)
    colnames(NucsLong) <- 1:dim(NucsLong)[[2]]
    NucsLong <- NucsLong[,colnames(Nucs)]
    saveRDS(NucsLong,NucsLongFN)
  }  
  }  
  
  expr.20 = quantile(SGD$Genes$expr, .2, na.rm  = TRUE )
  LongGenes = as.character(SGD$Genes$acc[which(width(SGD$Genes) > 2000 &
                                                 SGD$Genes$expr > expr.20)])
  MidGeneNucsId = which(mcols(NucRegions)$gene_pos  > 4 
                        & mcols(NucRegions)$acc %in% LongGenes)
  MidGenesNucs = ( colnames(Nucs) %in% MidGeneNucsId)
  
  print("Aggergating libraries")
  Nucs = AggregateNucs(Nucs)

  mod.nucs = ReadHistoneModeAtlas("Weiner-HisMod.csv", NucRegions, c() )
  m = mcols(mod.nucs)
  Inputs =  as.vector(exp(as.matrix(m["Input"])))
  NucInputs = as.vector(exp(as.matrix(m["Input"])))
  Inputs = Inputs[as.integer(colnames(Nucs))]
  names(Inputs) = colnames(Nucs)
  Inputs = Inputs/mean(Inputs)
  
  ColNames = colnames(Nucs)[Sample]
  NucRegions.select = resize(NucRegions[as.integer(ColNames)], width = 75, fix="center")
}

ComputeMarginalDistribution <- function( s, a1, tfactor = 1.1, threshold = .995) {
  N = Nuc(s,a1)
#  l1 = quantile(N, 1-threshold)*tfactor
  l1 = 0
  t1 = quantile(N, threshold)*tfactor
  pmax((Nuc(s,a1) -l1)/ (t1-l1), 0)
}

ComputeMarginalDistributionFactor <- function( s, a1, tfactor = 1.1, threshold = .995) {
  N = Nuc(s,a1)
  #  l1 = quantile(N, 1-threshold)*tfactor
  l1 = 0
  t1 = quantile(N, threshold)*tfactor
  1/(t1-l1)
}

PairwiseModelL <- function(N, err) {
#  sum(err/sqrt(N))
  -sum(N*log(N-err)-(N-err)- lfactorial(N))
}

PairwiseSlopeL <- function( N, Nmin, Nmax, a ) {
  Nmin = Nmin /a
  Nmax = Nmax /a
  Imax = N > 0 & Nmax > 0 & N > Nmax
  Imin = N > 0 & Nmin > 0 & N < Nmin
  PairwiseModelL( N[Imax], (N-Nmax)[Imax]) + PairwiseModelL( N[Imin], (N-Nmin)[Imin])
}

ComputePairwiseExpectations <- function (s, a1, a2, tfactor = 1.1, method = "L1",
                                         threshold = .995, Tmin = 0.05, symmetric = TRUE ) {
  N1 = ComputeMarginalDistribution(s, a1, tfactor, threshold )
  N2 = ComputeMarginalDistribution(s, a2, tfactor, threshold )
  N12 = Nuc(s,a1,a2)
  if( symmetric )
    N12 = N12 + Nuc(s,a2,a1)
#  l12 = quantile(N12, 1-threshold)*tfactor
  l12 = 0
  N12 = pmax(N12 - l12,0)
  r.exp = N1 * N2
  r.max = pmin(N1,N2)/N12
  #      hist(r.max,breaks="FD")
  r.min = pmax(N1+N2-1,0)/N12
  #      hist(r.min,breaks="FD")
  
  if( method == "percentile") {
    r.range = seq(0,1,by=0.005)
    r.max.q = quantile(r.max[r.exp>Tmin], r.range, na.rm=TRUE)
    r.min.q = quantile(r.min[r.exp>Tmin], 1-r.range, na.rm=TRUE)
    # find break tie between min and max violations
    #      plot(r.max.q,r.range, type="l")
    #      lines( r.min.q, r.range, col="red")
    if( any(r.min.q < r.max.q) ) {
      t12 = 1/r.max.q[which(r.min.q < r.max.q)[[1]]]
    } else {
      print(paste(s,a1,a2))
      print(r.min.q)
      print(r.max.q)
      t12 = 1/r.max.q[[20]]
    }
  } else {
    if( method == "ll") {
      ll <- function(x) { PairwiseSlopeL(N12, N12*r.max,N12*r.min, x) }
    } else
      if( method == "L1") {
        ll <- function(x) {
          Merr = N12-N12*r.max/x
          merr = N12*r.min/x - N12
          sum(Merr[Merr > 0 & !is.na(Merr)]) + sum(merr[merr > 0& !is.na(merr)])
        }
      } else
          abort(-1)

    opt = optimize(ll,lower=0.1*(1/max(N12)),upper=10)
#    print(opt)
    t12 = 1/opt$minimum
  } 
  list( l12 = l12, t12 = t12, 
        V.min = pmax(N1+N2-1,0), 
        V.exp = N1 * N2,
        V.max = pmin(N1,N2),
        V.obs = N12/t12
  )
}

ModelDistribution <- function( s, a1, a2 = "Input") {
  if( a2 == "Input" ) {
    D = ComputeMarginalDistribution(s,a1)
  } else {
    mod = ComputePairwiseExpectations(s,a1,a2, symmetric = TRUE)
    D = mod$V.obs
  }
  D
}

Plot2DProjectionGrid <- function( X, Y, 
                                  diagonal = TRUE, grid = 0.01, xmax = 1, legend = TRUE) {
  #    a = Project2DGrid(X,Y,grid)
  # a$N = log10(a$N)
  #    p = ggplot(a,aes(x=X,y=Y,color=N))
  #    p = p+geom_point(size=.5)
  #    p = p+scale_colour_gradient(low="blue",high="red")
  #    p
  I = X > 0.025 & Y > 0.025
  p = ggplot(data.frame(x  = X[I], y= Y[I]),aes(x=x,y=y))
  p = p+geom_bin2d(binwidth=grid)
  #    p = p+scale_fill_gradient(low="black",high="red")
  p = p+scale_fill_distiller(palette = "YlOrBr", direction=1, values=c(-0.1,.05,0.15,1))
  if( diagonal )
    p = p+geom_abline(slope=1,intercept=0,color="blue",size=.5)
  p = p + coord_equal(ratio=1)
  p = p+theme_classic()
  p = p + scale_x_continuous(expand=c(0,0),limits=c(0,xmax*1.1)) 
  p = p + scale_y_continuous(expand=c(0,0),limits=c(0,1.1)) 
  if( legend ) {
    p = p + theme(legend.title = element_text(size=10),legend.text = element_text(size=10))
  } else 
    p = p + theme(legend.position="none")
  p
}

library("cowplot")


CheckPlots <- function(s,a1,a2, Tmin = 0.05, method = "L1", threshold=.995,
                       title = paste0(s,"-",a1,"-",a2), Subset = NULL, legend = TRUE) {
  mod  = ComputePairwiseExpectations(s,a1,a2, method=method, Tmin = Tmin, threshold=threshold)
  N12 = mod$V.obs
  N1 = ComputeMarginalDistribution(s,a1,threshold=threshold)
  N2 = ComputeMarginalDistribution(s,a2,threshold=threshold)
  p1 = Plot2DProjectionGrid(N1,N12, legend=legend)
  p1 = p1 + labs(x = a1, y = paste0(a1,"-",a2))
  p2 = Plot2DProjectionGrid(N2,N12, legend=legend)
  p2 = p2 + labs(x = a2, y = paste0(a1,"-",a2))
  p3 = Plot2DProjectionGrid((N1+N2),N12, diagonal = FALSE, xmax = 2, legend=legend)
  p3 = p3 + labs(x = paste0(a1,"+",a2), y = paste0(a1,"-",a2))
  p3 = p3 + geom_abline(slope = 1, intercept = -1,color="blue", size=.5)
  p4 = Plot2DProjectionGrid(mod$V.exp,N12, legend=legend)
  p4 = p4 + labs(x = paste0(a1,"*",a2), y = paste0(a1,"-",a2))
  if( is.null(Subset) ) {
    plot_grid(p1,p2,p3,p4, ncol=1)
  } else {
    p5 = Plot2DProjectionGrid(N1[Subset],N12[Subset], legend=legend)
    p5 = p5 + labs(x = a1, y = paste0(a1,"-",a2))
    p6 = Plot2DProjectionGrid(N2[Subset],N12[Subset], legend=legend)
    p6 = p6 + labs(x = a2, y = paste0(a1,"-",a2))
    p7 = Plot2DProjectionGrid((N1[Subset]+N2[Subset]),N12[Subset], diagonal = FALSE, xmax = 2, legend=legend)
    p7 = p7 + labs(x = paste0(a1,"+",a2), y = paste0(a1,"-",a2))
    p7 = p7 + geom_abline(slope = 1, intercept = -1,color="blue", size=.5)
    p8 = Plot2DProjectionGrid(mod$V.exp[Subset],N12[Subset], legend=legend)
    p8 = p8 + labs(x = paste0(a1,"*",a2), y = paste0(a1,"-",a2))
    plot_grid(p5,p1,p6,p2,p7,p3,p8,p4,ncol = 2)
  }
}

