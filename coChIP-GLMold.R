##
## To describe a GLM problem we need:
##
## Length of list of unknown per site
## Matrix mapping from unknowns to prediction in each file
##

library(optimx)
library(lbfgsb3)

##
## mutli-design
##
## A = A_1,...,A_n s.t. Y_i = theta^t A_i x
##
ccMultiDesign2RowDesign <- function(multiDesign, par, Nus) {
  do.call(rbind, mapply(function(A,nu) nu*(t(par)%*%A),multiDesign, Nus, SIMPLIFY = FALSE) )
}

ccMultiDesign2ColDesign <- function(multiDesign, Xs, Nus) {
  do.call(rbind, mapply(function(A,nu) nu*(Xs %*% t(A)), multiDesign, Nus, SIMPLIFY = FALSE))
}

##
## Approx of hinge function
##
hinge.rad = .1
#hinge.rad2 = hinge.rad**2
#hinge.xc = 1 - hinge.rad*tan(pi/8)
#hinge.tresh = hinge.xc +hinge.rad/sqrt(2)

ccGLMhinge <- function(x) {
   y = x
   y[y<0] = 0
#   y[x > hinge.xc & x < hinge.tresh] = hinge.rad - sqrt(hinge.rad2 - (x[x > hinge.xc & x < hinge.tresh] - hinge.xc)**2)
   y[x > -hinge.rad & x < +hinge.rad] = (x[x > -hinge.rad & x < +hinge.rad]+hinge.rad)**2/(4*hinge.rad)
   return(y)
}

ccGLMhingeGrad <- function(x) {
  y = 1*(x>0)
#  y[x > hinge.xc & x < hinge.tresh] = (x[x > hinge.xc & x < hinge.tresh] - hinge.xc)/sqrt(hinge.rad2 - (x[x > hinge.xc & x < hinge.tresh] - hinge.xc)**2)
  y[x > -hinge.rad & x < hinge.rad] = (x[x > -hinge.rad & x < +hinge.rad]+hinge.rad)/(2*hinge.rad)
  return(y)
}



##
## Do row estimate
##
ccGLMSingleRowFit.glm <- function( Y, Design, initX ) {
  g = glm(Y ~ Design - 1, family = poisson(link="identity"), start = initX )
  x = coef(g)
  x[is.na(x)] = 0
  return(c(x, deviance(g)))  
}

ccGLMSingleRowFit <- function( Y, Design, initX, penalty = 10000 ) {
  neg.log.like <- function (x) {
#    print("ll")
#    print(x)
    z = Design %*% x
    l = -sum( dpois(Y, z, log=TRUE) )
#    l = l + penalty*(1-sum(x))**2
    sx = sum(x)
#    if( sx > 1 )
#      l = l + penalty*(sx-1)**2
    l = l +penalty*ccGLMhinge(sx-1)
#    print(paste( do.call(paste,as.list(x)), l))
    l
  }
  grad.log.like <- function( x ) {
#    print("grad")
#    print(x)
    z = Design %*% x
    dl = -(Y / t(z) - 1)%*% Design
    #dl = dl - 2*x*penalty*(1-sum(x))
    sx = sum(x)
#    if(sx > 1)
#      dl = dl + 2*penalty*(sx -1)
    dl = dl + penalty*ccGLMhingeGrad(sx-1)
    return(dl)
  }
#  g = optimx(initX, neg.log.like, gr=grad.log.like, lower = 0.00001, method="L-BFGS-B")
#  return(list(as.vector(coef(g)), g$value))
  g = lbfgsb3(initX, neg.log.like,grad.log.like,lower = 0.00001,control=list(iprint=-1))
#  print( (g))
  return(list(as.vector(g$prm), g$f))
}

ccGLMGeneralLogLikelihood <- function( counts, multiDesign, par, Nus, Xs, penalty = 1000 ) {
  print("ll")
  Design = ccMultiDesign2RowDesign(multiDesign, par, Nus)
  z = Design %*% t(Xs)
  if( min(z) < 0 )
    return(Inf)
  sx = rowSums(Xs)
  l = -sum(dpois(counts,z,log=TRUE)) + penalty*sum(ccGLMhinge(sx - 1))
  return(l)
}

ccGLMGeneralLogLikelihoodGrad <- function( counts, multiDesign, par, Nus, Xs, penalty = 1000 ) {
  print("gl")
  Design = ccMultiDesign2RowDesign(multiDesign, par, Nus)
  z = Design %*% t(Xs)
  sx = rowSums(Xs)
  m = dim(counts)[[2]]
  k = dim(Design)[[2]]
  dl = -t(counts /z -1) %*% Design + 
        matrix(penalty*ccGLMhinge(sx - 1),nr=m,nc=k)
  return(dl)
}

ccGLMFitRowsCombined <- function( counts, multiDesign, par, Nus, prevXs = NULL,   penalty = 1000 ) {
  Design = ccMultiDesign2RowDesign(multiDesign, par, Nus)
  m = dim(counts)[[2]]
  k = dim(Design)[[2]]
  
  ll <- function( p ) {
    ccGLMGeneralLogLikelihood(counts, multiDesign, par, Nus, matrix(p, nr = m, nc = k), penalty )
  }  
  lg <- function( p ) {
    ccGLMGeneralLogLikelihoodGrad(counts, multiDesign, par, Nus, matrix(p, nr = m, nc = k), penalty )
  }  

  if( is.null( prevXs) )
    initXs = matrix(1/k,nc = k, nr = m)
  else
    initXs = prevXs
  
  g = lbfgs(ll,lg,as.vector(initXs))
  print(g)
  matrix(g$value, nc=k, nr=m)
}

##
## Fit region based unknown
##
ccGLMFitRows <- function( counts, multiDesign, par, Nus, prevXs = NULL ) {
  print("Fit rows")
  Design = ccMultiDesign2RowDesign(multiDesign, par, Nus)
  print(Design)
  m = dim(counts)[[2]]
  k = dim(Design)[[2]]
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
    r = ccGLMSingleRowFit(Y, Design, start)
    x = r[[1]]
    d = d+r[[2]]
#      },
#      error = function( err ) {
#        print(paste("Error ",err, "i = ", i, "Y = ", do.call(paste,as.list(Y))))
#      })
    Xs[i,] = x
  }
  print(d)
  # renormalize 
  #return(Xs/max(Xs))
  Xs
}

ccGLMFitCols <- function( counts, multiDesign, Xs, Nus) {
  print("Fit columns")
  Design = ccMultiDesign2ColDesign(counts, multiDesign, Xs, Nus)
  Y = as.vector(counts)
  g = glm(Y ~ Design -1, family = poisson(link="identity"))
  print(deviance(g))
  par = coef(g)
}

ccGLMFitNus <- function( counts, multiDesign, par, Xs ) {
  # First normalize Xs to be close to probability
  m = median(rowSums(Xs))
  Xs = Xs/m
  
  # Compute counts vs predictions (without Nu) in each experiment
  Cs = unlist(lapply(1:length(multiDesign), function(i) sum(counts[i,],na.rm = TRUE) ))
  Ps = unlist(lapply(1:length(multiDesign), 
               function(i)  sum( Xs%*%t(par %*% multiDesign[[i]]))))
  print( paste( "Counts:", do.call(paste,as.list(Cs))))
  print( paste( "Pred:", do.call(paste,as.list(Ps))))
  
  # Set Nus to make the prediction match the sum
  lapply(1:length(multiDesign), function(i) Cs[i]/Ps[i])
}

## counts - matrix of counts (region/nucs/etc). Rows - samples, Col - regions
## metaDesign - n matrixs so that Y_i = par^t A_i X

ccGLMFit <-function( counts, metaDesign, initialNuInflation = 10, 
                     initialPar = NULL, 
                     updatePar = FALSE,
                     initialNu = NULL,
                     iter = 5) {
  n = dim(counts)[[1]]  # number of samples
  m = dim(counts)[[2]]  # number of regions
  k = dim(metaDesign[[1]])[[1]] # number of parameters
  l = dim(metaDesign[[1]])[[2]] # number of unknowns
  
  if( n != length(metaDesign) ) {
    print(paste("Dimension mismatch ",n,"vs",length(metaDesign)))
    return(0)
  }
  
  ## Initial guess for count sample size
  if( is.null(initialNu) )
    Nus <- rowSums(counts) * initialNuInflation
  else
    Nus <- initialNu
  
  if( is.null(initialPar) ) {
    if( !updatePar )
      print("Can't have no initial parameters and no update at the same time");
    initialPar = runif(k)
  }
  par = initialPar
  for( i  in 1:iter ) {
    Xs = ccGLMFitRows(counts, metaDesign, par, Nus)
    Nus = ccGLMFitNus(counts, metaDesign, par, Xs)
    if( updatePar ) {
      par = ccGLMFitCols(counts, metaDesign, Xs, Nus)    
      Nus = ccGLMFitNus(counts, metaDesign, par, Xs)
      print(par)
    }
    print(Nus)
  }
  list( par, Xs, Nus )  
}