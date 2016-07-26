## Old and inefficent 

AvgRegionsSimple <- function( cov, rs, width=1500 ) {
  r.by.chrom = split(ranges(rs), seqnames(rs))
  sv = replicate(width,0)
  m = 0
  sumv <- function(v) {
    #    sv <<- rowSums(cbind(as.numeric(v),sv))
    #sv <<- mapply(sum,as.numeric(v),sv)
    sv <<- as.numeric(v)+sv
    1
  }
  for( x in labels(r.by.chrom)) {
    m  = length(cov[[x]])
    I = start(r.by.chrom[x]) > 0 & end(r.by.chrom[x]) < m
    vs = Views(cov[x], r.by.chrom[x][I])
    lapply(vs[[1]],sumv)
  }
  n = length(ranges(rs))
  return( sv )
}

AvgRegionsTwoStrands <- function(cov, rs, plus, width=1500) {
  ap = AvgRegionsSimple(cov,rs[plus],width)
  an = AvgRegionsSimple(cov,rs[!plus],width)
  a = ap + an[length(an):1]
  return(a/length(plus))
}

AvgRegionsSubgroup <- function(cov, rs, plus, A, width=1500)
{
  I = is.element(rs$acc, A)
  AvgRegionsTwoStrands( cov, rs[I], plus[I], width )
}
AvgRegionsSubgroups <- function(cov, rs, plus, subs, width=1500)
{
  A = lapply(subs, function(A) AvgRegionsSubgroup(cov,rs,plus,A,width) )
  names(A) = names(subs)
  return(A)
}
