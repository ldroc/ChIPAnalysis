FragmentBoundaries = c( 100, 200, 250, 300, 350, 400, 450, 500, 600, 700, 1000 )

CoverageByLength <- function( GR, Lengths ) {
  n = length( Lengths )
  W = width(GR)
  
  SelectRange <- function (i) {
    I = (W < Lengths[i+1]) & (W >= Lengths[i]); 
    return(coverage(GR[I]))
  }
  l = lapply( 1:(n-1), SelectRange ) 
  names(l) = unlist(lapply(1:(n-1), function(i) paste(Lengths[i], "-", Lengths[i+1])))
  return(l)
}

CollectRegionsBySize <- function( LenCov, regions, width = 1500) {
  lapply(LenCov, function(c) CollectRegions(c, regions, width) )
}

AvgRegionsLenBySubgroups <- function( Regs, rs, subs ) {
  lapply(Regs, function(c) AvgRegionsSubgroups(c,rs,subs))
}

PlotSizeCoverage <- function (AA, label="TSS", pos=500, file = NULL, path=NULL ) {
  if( !is.null(file) ) {
    if( !is.null(path) )
      out = paste0(path,"/", BaseFileName(file),"-Meta.pdf")
    else
      out = paste0(file_path_sans_ext(file),"-Meta.pdf")
    pdf(out,paper="a4",pointsize = 10)
  }
  
  # there must be a better way to do this...
  m = max(unlist(lapply(AA, function(x) max(unlist(lapply(x,max))))))
  ylim = c(0,m)
  ns = rev(names(AA[[1]]))
  print(ns)
  par(mfrow = c(length(ns),1))
  parparam =  par(mar=c(3.1,4.1,1.5,2.1))
  
  lapply(ns, 
         function(l)  PlotCovergeGroups(lapply(AA, function(a) a[[l]]), 
                                        main=l, ylim = ylim,
                                        label=label, pos=pos) )
  par(parparam,mfrow=c(1,1))
  if( !is.null(file) )
    dev.off()
    
}

PlotMultiSizeCoverage <- function( Files, AA, path = NULL, pos=500, label="TSS") {
  lapply(1:length(Files), function(i) PlotSizeCoverage( AA[[i]], file = Files[i], 
                                                      path=path, pos=pos, label=label))
}
