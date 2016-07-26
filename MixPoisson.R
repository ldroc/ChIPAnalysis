FitNonZeroNB <- function( ObsStat ) {

#  print(ObsStat)
  L = length(ObsStat)
  M = as.integer(names(ObsStat))
  N = sum(ObsStat)
  
  ll <- function( x )
  {
    s = x[[1]]
    p = x[[2]]
    -sum(ObsStat * dnbinom(M, s, p, log = TRUE )) - N*log(1-dnbinom(0, size=s, prob=p, log=FALSE))
  }

  x = c(1, .5)
  tryCatch( {
    opt = optim(x, ll )
    return(opt$par)
    
  }, error = function(e) print(e) )
  return(x)
}


FitNonZeroPoisson <- function( ObsStat ) {
  N = sum(ObsStat)
  M = max(as.integer(names(ObsStat)))
  E = sum(ObsStat * as.integer(names(ObsStat))) / N
  ll <- function( p ) {
    -(E * log(p) - p - log(1-exp(-p)))
  }
  tryCatch( {
    opt = optimize(f = ll, interval = c(0,M))
    return(opt$minimum)
  }, error = function(e) print(e) )
  return(1)
}

MixPoisson <- function(x,lambda,p) {
  y = x*0
  for( i in 1:length(p) ){
    y = y+ p[i]*dpois(x,lambda*i)
  }
  return(y)
}

MixPoissonAvg <- function(lambda,p) {
  
  y = 0
  for( i in 1:length(p) ){
    y = y+ p[i]*lambda*i
  }
  return(y)
}

MixPoissonLL <- function(Obs, lambda, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12) {
  x = 1:length(Obs)
  p = c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
  z = sum(exp(p))
  q = MixPoisson( x, lambda, exp(p)/z )
  t = 1-MixPoisson( 0, lambda, exp(p)/z )
#  print(Obs)
  qll = -sum( Obs * log(q/t) )
}

FitMixPoisson <- function( ObsStat )
{
  i = list( lambda = FitNonZeroPoisson(ObsStat)/2, p1=10, p2=0, p3=0, p4=0,p5=1,p6=0,p7=0,p8=0,p9=0,p10=0,p11=0,p12=0)
#  ll <- function( lambda, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12) {
#  return(MixPoissonLL(ObsStat, lambda, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 ))
#}
  #  est <- mle(minuslogl=ll, start=i )
  #estparam = coef(est)
  
    ll <- function( p ) {
    
      return(MixPoissonLL(ObsStat, p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13]  ))
  }

#  est <- optim(i,ll,hessian = FALSE,method="SANN", control=list(maxit=10000))
    estparam = i
    tryCatch( {
      est <- optim(i,ll,control=list(maxit=10000))
      #      print(est$value)
      estparam = est$par
      p = c(estparam["p1"], estparam["p2"], estparam["p3"], estparam["p4"], 
            estparam["p5"], estparam["p6"], estparam["p7"],estparam["p8"],
            estparam["p9"], estparam["p10"], estparam["p11"], estparam["p12"])
      p = exp(p)
      p = p/sum(p)
      lambda = estparam["lambda"]
      return(list( lambda = lambda, p = p))
      
    }, error = function(e) print(e) )
#  print(lambda)
#  plot(table(dups),type='h'); plot(1:31, MixPoisson(1:31,lambda, p), type = 'h')
    return(list(lambda = 1, p = c(1,0,0,0)))
}