##
##
##

##
## Func
##
## Design(params) value: matrix[numPar,numX]
## DesignGrad(params) value: list of numPar matrices, each [numPar,numX]
## Init(counts) value: list(Xs = ,par = )
## numX value: number of unknowns
## numPar value: number of parameters
##
library(optimx)
library(lbfgsb3)

Xeps = 10**-6

##
## Approx of hinge function
##
hinge.rad = .1

ccGLMhinge <- function(x) {
   y = x
   y[y<0] = 0
   y[x > -hinge.rad & x < +hinge.rad] = (x[x > -hinge.rad & x < +hinge.rad]+hinge.rad)**2/(4*hinge.rad)
   return(y)
}

ccGLMhingeGrad <- function(x) {
  y = 1*(x>0)
  y[x > -hinge.rad & x < hinge.rad] = (x[x > -hinge.rad & x < +hinge.rad]+hinge.rad)/(2*hinge.rad)
  return(y)
}

##
## Example of parameter function
##
ccPairIPDesign <- function( params ) {
  a.spe = params[[1]]
  b.spe = params[[2]]
  a.bg = params[[3]]
  b.bg = params[[4]]
  D = matrix(c( a.spe, a.spe, a.bg, a.bg,  # a-input
                b.spe, b.bg, b.spe, b.bg,   # b-input
                a.spe*a.spe, a.spe*a.spe, a.bg*a.bg,  a.bg*a.bg, # a-a
                b.spe*b.spe, b.bg*b.bg, b.spe*b.spe, b.bg*b.bg, # b-b
                a.spe*b.spe, a.spe*b.bg, a.bg*b.spe, a.bg*b.bg, # a-b
                a.spe*b.spe, a.spe*b.bg, a.bg*b.spe, a.bg*b.bg),  # b-a
             nr = 6, nc = 4, byrow=TRUE)
  return(D)
} 

ccPairIPDesignGrad <- function( params ) {
  a.spe = params[[1]]
  b.spe = params[[2]]
  a.bg = params[[3]]
  b.bg = params[[4]]
  Grad = list( 
    matrix(c( 1, 1, 0, 0,  # a-input
              0, 0, 0, 0,   # b-input
              a.spe, a.spe, 0,  0, # a-a
              0, 0, 0, 0, # b-b
              b.spe, b.bg, 0, 0, # a-b
              b.spe, b.bg, 0, 0), # b-a
           nr = 6, nc = 4, byrow=TRUE),
    matrix(c( 0, 0, 0, 0,  # a-input
              1, 0, 1, 0,   # b-input
              0, 0, 0, 0, # a-a
              b.spe, 0, b.spe, 0, # b-b
              a.spe, 0, a.bg, 0, # a-b
              a.spe, 0, a.bg, 0),  # b-a
           nr = 6, nc = 4, byrow=TRUE),
    matrix(c( 0, 0, 1, 1,  # a-input
              0, 0, 0, 0,   # b-input
              0, 0, a.bg,  a.bg, # a-a
              0, 0, 0, 0, # b-b
              0, 0, b.spe, b.bg, # a-b
              0, 0, b.spe, b.bg),  # b-a
           nr = 6, nc = 4, byrow=TRUE),
    matrix(c( 0,0,0,0,  # a-input
              0, 1, 0, 1,   # b-input
              0, 0, 0, 0,  # a-a
              0, b.bg, 0, b.bg, # b-b
              0, a.spe, 0, a.bg, # a-b
              0, a.spe, 0, a.bg),  # b-a
           nr = 6, nc = 4, byrow=TRUE) )
  return(Grad)
} 

ccPairIPInit <- function( counts, input = NULL, spec = .5, bg = 0.01 ) {
  threshold = 0.995
  ## Initialize using the input samples and assume there is joint activity
  n = dim(counts)[[2]]
  Xs = matrix(0, nr = n, nc = 4)

  if( is.null(input) )
    input = rep(1, n)
  
  I = input / quantile(input,threshold, na.rm=TRUE)
  A = counts[1,] / I
  A = A / quantile(A,threshold, na.rm=TRUE)
  B = counts[2,] / I
  B = B / quantile(B,threshold, na.rm=TRUE)

  A = pmin(A,1)
  B = pmin(B,1)
  Xs[,1] = A * B
  Xs[,2] = A * (1-B)
  Xs[,3] = (1-A) * B
  Xs[,4] = (1 - rowSums(Xs[,1:3]))
  Xs = Xs * I
  
  par = c( spec, spec, bg, bg)
  
  return( list(Xs = Xs, par = par))
}


ccGLMGeneralDesignMat <- function( Func, Nus) {
  Design = Func$Design(par) 
  d1 = dim(Design)[[1]]
  d2 = dim(Design)[[2]]
  Design = Design * matrix( Nus, nr = d1, nc = d2)
  return(Design)
}

ccGLMGeneralLogLikelihood <- function( counts, Func, par, Nus, Xs, penalty = 1000 ) {
  print("ll")
  Design = ccGLMGeneralDesignMat( Func, Nus )
  z = Design %*% t(Xs)
  sx = rowSums(Xs)
  l = -sum(dpois(counts,z,log=TRUE)) + penalty*sum(ccGLMhinge(sx - 1))
  return(l)
}

ccGLMGeneralLogLikelihoodGradXs <- function( counts, Func, par, Nus, Xs, penalty = 1000 ) {
  print("gl-x")
  Design = Func$Design(par) 
  d1 = dim(Design)[[1]]
  d2 = dim(Design)[[2]]
  Design = Design * matrix( Nus, nr = d1, nc = d2)
  z = Design %*% t(Xs)
  sx = rowSums(Xs)
  m = dim(counts)[[2]]
  k = dim(Design)[[2]]
  dl = -t(counts /z -1) %*% Design + 
        matrix(penalty*ccGLMhinge(sx - 1),nr=m,nc=k)
  return(dl)
}

ccGLMGeneralLogLikelihoodGradParNu <- function( counts, Func, par, Nus, Xs, penalty = 1000 ) {
  print("gl-p-n")
  Design = Func$Design(par) 
  d1 = dim(Design)[[1]]
  d2 = dim(Design)[[2]]
  NuMat = matrix( Nus, nr = d1, nc = d2)
  z = Design %*% t(Xs)
  z = z * NuMat
  NuGrad = (colSums(counts) - colSums(z))/Nus
  
  m = dim(counts)[[2]]
  k = dim(Design)[[2]]

  M = counts/z -1
  GG = Func$Grad(par)
  parGrad = lapply(GG, function(G) sum(M*(G %*% t(Xs))))
  
  return( c(parGrad, NuGrad) )
}

ccGLMXsFit <- function( Y, Func, par, initXs, penalty = 10000 ) {
  m = dim(Xs)[[1]]  
  n = dim(Xs)[[2]]
  neg.log.like <- function (x) {
    ccGLMGeneralLogLikelihood( Y, Func, par, Nus, matrix(x, nr=m, nc = n), penalty ) 
  }
  grad.log.like <- function( x ) {
    as.vector(ccGLMGeneralLogLikelihoodGradXs( Y, Func, par, Nus, matrix(x, nr=m, nc = n), penalty ))
  }
  g = lbfgsb3(as.vector(initXs), neg.log.like,grad.log.like,lower = Xeps, control=list(iprint=-1))
  #  print( (g))
  return(list(Xs = matrix(g$prm, nr=m, nc=n), ll = g$f))
}

##
## Fit region based unknown
##
ccGLMFitRows <- function( counts, Func, par, Nus, prevXs = NULL ) {
  print("Fit rows")
  m = dim(counts)[[2]]
  k = Func$numX
  Xs = matrix(nrow=m,ncol=k)
  d = 0
  # run over regions
  for( i in 1:m )  {
    Y = counts[,i]
    x = rep(0,k)
#    tryCatch( {
    start = rep(1/k,length=k)
    if( !is.null( prevXs ))
      start = prevXs[i,]
    r = ccGLMXsFit(Y, Func, par, Nus, start)
    x = r[[1]]
    d = d+r[[2]]
    Xs[i,] = x
  }
  print(d)
  return(Xs)
}

ccGLMFitParNus <- function( counts, Func, par, Nus, Xs) {
  print("Fit par Nus")
  l = 1:length(par)
  m = (1:length(Nus))+length(par)
  
  neg.log.like <- function (x) {
    ccGLMGeneralLogLikelihood( Y, Func, x[l], x[m], Xs, penalty ) 
  }
  grad.log.like <- function( x ) {
    ccGLMGeneralLogLikelihoodGradParNu( Y, Func,  x[l], x[m], Xs, penalty )
  }
  g = lbfgsb3(c(par,Nus), neg.log.like,grad.log.like,lower = 0, control=list(iprint=-1))
  
  return(list(par = g$prm[l],Nus = g$prm[m], ll = g$f))
}

## counts - matrix of counts (region/nucs/etc). Rows - samples, Col - regions
## metaDesign - n matrixs so that Y_i = par^t A_i X

ccGLMFit <-function( counts, Func, iter = 5 ) {
  n = dim(counts)[[1]]  # number of samples
  m = dim(counts)[[2]]  # number of regions
  k = Func$numPar # number of parameters
  l = Func$numX # number of unknowns
  
  init <- Func$Init(counts)
  par = init$par
  Xs = init$Xs

  for( i  in 1:iter ) {
    Xs = ccGLMFitRows(counts, Func, par, Nus, Xs)
    pn = ccGLMFitParNus(counts, Func, par, Nus, Xs)
    par = pn$par
    Nus = pn$Nus
    print(paste("Iteration", i, " ll = ", pn$ll))
    print(par)
    print(Nus)
  }
  list( par = par, Xs = Xs, Nus = Nus )  
}