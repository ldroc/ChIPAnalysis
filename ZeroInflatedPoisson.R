ZIP.sample <- function( p, lam, Ns ) {
  rpois(length(Ns), Ns*lam)*rbinom(length(Ns),1,p)
}

ZIP.pdf <- function( p, lam, cs, Ns) {
  prob <- p*((lam*Ns)**cs)*exp(-lam*Ns)/factorial(cs)
  prob[cs==0] = prob[cs==0]+1-p
  prob
}

ZIP.loglikehood <- function( p, lam, cs, Ns) {
  L = log(p) + cs*log(lam) + cs*log(Ns) -lam*Ns -lfactorial(cs)
  Z = cs == 0
  L[Z] = log(exp(L[Z]) + 1-p)
  sum(L)
}


ZIP.loglikehood.gradient <- function( p, lam, cs, Ns) {
  L = log(p) + cs * log(lam) + cs*log(Ns) -  lam*Ns - lfactorial(cs)
  Ll = cs/lam - Ns
  Z = cs==0
  Ll[Z] = Ll[Z]*exp(L[Z])/(exp(L[Z]) + 1-p)
  Lp = as.vector(rep(1/p,length(cs)))
  Lp[Z] = (Lp[Z]*exp(L[Z])-1)/(exp(L[Z]) + 1-p) 
  c(sum(Lp), sum(Ll))
}

ZIP.estimate <- function( cs, Ns ) {
  p = sum( cs > 0 ) / length(cs)
  lam = mean(cs[cs>0])/mean(Ns[cs>0])
  
  eps = 1e-7
  opt = optim(par = c(p = p,lam = lam), 
        function(pvec) {-ZIP.loglikehood(pvec[1], pvec[2], cs, Ns)},
        function(pvec) {-ZIP.loglikehood.gradient(pvec[1], pvec[2], cs, Ns)},
        lower=c(eps,eps),
        upper=c(1-eps,Inf)
  )
  opt$par
}

ZIP.testdata <- function(p,lam,n, mu) {
  Ns = rpois(n,mu)
  cs = ZIP.sample(p,lam,Ns)
  list( cs = cs, Ns = Ns)
}

p = 0.3
lam = 0.0005
mu = 10000
n = 1000
testdata <- ZIP.testdata(p,lam,n,mu)
ggplot( data.frame(x = testdata$cs), aes(x=x)) + geom_histogram(binwidth=1)

ZIP.estimate(testdata$cs, testdata$Ns)
