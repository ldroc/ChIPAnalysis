QCGRanges <- function( BAMFile, dat, output = NULL ) {
  
  BAMFile = BaseFileName(BAMFile)
  dat.uniq = unique(dat)
  
  if(length(dat.uniq) < 100)
    return(0)
  
  if( !is.null(output) )
    png(output, height = 1024, width=1024)
  
  par(mfrow = c(2,1))
  
  # Build a histogram of fragment lengths
  title = paste0('Fragment lengths (',BAMFile, ")  ", length(dat.uniq), " unique frags")
  hist( width(dat.uniq), main = title , 
        xlab = 'fragement length', 
        xlim = c(0,1000),
        breaks="FD")
  
  # count duplicate segments
  dups <- countOverlaps(dat.uniq, dat, type = "equal")
#  dup.hist = rle(sort(dups))
  
  m = max(dups)
  dupStat = table(dups)
  p = FitMixPoisson(dupStat)
  N = sum(length(dat.uniq))
  lambda = as.integer(p$lambda*100)/100.0
  e = sum(unlist(lapply(1:length(p$p), function (i) i*p$p[[i]])))
  e = as.integer(e*100)/100.0
  pge0 = 1-MixPoisson(0,p$lambda, p$p)
  duptitle = paste0('Fragment duplicates (',BAMFile, ") ",length(dat), 
                   " frags.\n Lambda = ", lambda, ", Est N = ", as.integer(N/pge0) )
  hist( dups, main = duptitle, xlab = 'numbers', breaks=((-0.5:(m+.5))),
        xlim = c(0,max(dups)))
  
  pge0 = 1-MixPoisson(0,p$lambda, p$p)
  lines(-.55:(m-.45),MixPoisson(0:m,p$lambda, p$p)*N/pge0,col="red",type="s")
  l = FitNonZeroPoisson(dupStat)
  pge0 = 1 - dpois(0,l)
  lines(-.6:(m-.4),dpois(0:m,l)*N/pge0,col="blue",type="s")
#  ll = FitNonZeroNB(dupStat)
#  pge0 = 1 - dnbinom(0,ll[[1]], ll[[2]])
#  lines(-.57:(m-.43),dnbinom(0:m,ll[[1]], ll[[2]])*N/pge0,col="purple",type="s")
  
  
  #  sim = rle(sort(rnbinom(N, 
  #                         size = p$estimate["size"], 
  #                        mu = p$estimate["mu"])))
  #  lines(sim$values, sim$lengths, col="red")
  #sim = rle(sort(rpois(N, lambda = p$estimate["mu"])))
  #lines(sim$values, sim$lengths, col="blue")
  
  if( !is.null(output) )
    dev.off()
  par(mfrow = c(1,1))
}

QCSummaryRanges <- function(dat, SizeRanges = c(1,1000) ) {
  tmp <- function( dat ) {
    N = length(dat)
    dat.uniq = unique(dat)
    dups <- countOverlaps(dat.uniq, dat, type = "equal")
    dupStat = table(dups)
    Nuniq = length(dat.uniq)
    lambda1 = FitNonZeroPoisson(dupStat)
    Total1 = Nuniq /(1 - exp(-lambda1))
    p = FitMixPoisson(dupStat)
    lambda = p$lambda * sum(unlist(lapply(1:length(p$p), function (i) i*p$p[[i]])))
    Total = Nuniq /(1 - MixPoisson(0, p$lambda, p$p))
    list( N = N, Nuniq = Nuniq, Dup = N/Nuniq, EstTotal = Total1, EstTotalMix = Total, Lambda = lambda1, MixLambda = lambda )
  }
  
#  print(lapply(1:(length(SizeRanges)-1), function(i) c(SizeRanges[[i]],SizeRanges[[i+1]])))
  Q = lapply(1:(length(SizeRanges)-1), function(i) tmp(dat[width(dat) >= SizeRanges[[i]] & width(dat) < SizeRanges[[i+1]]]))
  names(Q) = lapply(1:(length(SizeRanges)-1), function(i) paste0(SizeRanges[[i]],"-",SizeRanges[[i+1]]))
  Q
}

QCGRangesMultileFiles <- function( rs, path=NULL ) {
  tf <- function( r ) {
    f = r$BAMFile
    if( !is.null(path) )
      out = paste0(path,"/", BaseFileName(r$BAMFile),"-QC.png")
    else
      out = paste0(file_path_sans_ext(r$BAMFile),"-QC.png")
    QCGRanges( f, r$GR, output = out )
  }
  sapply(rs, tf)
}


BaseFileName <- function( fname ) {
  x = file_path_sans_ext(fname)
  y = strsplit(x,"/")[[1]]
  n = length(y)
  z = y[n]
  return(z)
}
