

Expr = SGD$Genes$expr
names(Expr) = SGD$Genes$name
Expr = Expr[!is.na(Expr) & Expr > 0]
Expr = Expr/sum(Expr)
N.genes = length(Expr)

SampleTranscriptom <- function( N ) {
  rpois(N.genes, Expr*N)
}

SampleCapture <- function(Trans, p) {
  rbinom(N.genes, Trans, p)
}

SampleSeq <- function( C, lambda) {
  rpois(N.genes, C*lambda)
}

CreateCaptureGraph <- function() {
  df = data.frame(x = c(0), y=c(0), lambda=c(0))
  ggplot(df,aes(x=lambda,y=y))
}

AddCaptureGraph <- function(plot,p, lambda, col, N = 15000, n =100) {
  A = sapply(lambda, function(l) replicate(n,{S = SampleSeq(SampleCapture(SampleTranscriptom(N),p),l); c(sum(S), sum(S>0))} ), simplify = "array")
  B = do.call(cbind,lapply(1:dim(A)[3], function(i) A[,,i]))
  df = data.frame(x = B[1,],y=B[2,],lambda=rep(lambda, each=n))
  plot+geom_point(data=df, shape=16,alpha = 0.05, color=col)+geom_smooth(data=df,se=TRUE,color=col)
}

if(0) {
  p = CreateCaptureGraph()
  probs = c(.1,.2,.4,.8)
  for( i in 1:length(probs) )
    p = AddCaptureGraph(p,probs[i],seq(.1,4,by=.1),col=i)
}

if(0) {
  Ps = seq(0,.5,by=0.01)
  df = data.frame(p = c(0), zero = c(0),n=c(0))
  p = ggplot(df, aes(x = p, y = zero,color=n))
  for( i in 0:6 ) {
    n = 2**i
#  for( n in seq(1,8,by=1) ) {
    df = data.frame(p = Ps, zero = (1-Ps)**n, n = rep(format(n,width=3),length(Ps)))
    p = p + geom_line(data=df,size=2)
  }
  p = p + scale_color_brewer(palette = "Oranges")
  p = p+theme_minimal()
  p = p + labs(title="Prob of Zero given capture", x="Capture Efficiency", y = "Prob(zero)")
}


if(0) {
  Ps = seq(0,.5,by=0.01)
  df = data.frame(p = c(0), zero = c(0),n=c(0))
  p = ggplot(df, aes(x = p, y = zero,color=n))
  for( i in 0:6 ) {
    n = 2**i
    #  for( n in seq(1,8,by=1) ) {
    df = data.frame(p = Ps, zero = Ps*n, n = rep(format(n,width=3),length(Ps)))
    p = p + geom_line(data=df,size=2)
  }
  p = p + scale_color_brewer(palette = "Oranges")
  p = p+theme_minimal()
  p = p + labs(title="Expectation given capture", x="Capture Efficiency", y = "Expected count")
}

if(1) {
  N = 300000
  Ps = seq(0,.5,by=0.01)
  df = data.frame(p = c(0), zero = c(0),n=c(0))
  p = ggplot(df, aes(x = p, y = zero,color=n))
  for( i in 0:6 ) {
    n = 2**i
    #  for( n in seq(1,8,by=1) ) {
    df = data.frame(p = Ps, zero = sapply(Ps, function(prob) 10**6*median(rbinom(10000,n,prob)/rbinom(1000,N,prob))), n = rep(format(n,width=3),length(Ps)))
    p = p + geom_line(data=df,size=2)
  }
  p = p + scale_color_brewer(palette = "Oranges")
  p = p+theme_minimal()
  p = p + labs(title="Median TPM given capture", x="Capture Efficiency", y = "TPM")
  p.medianTPM = p
  df = data.frame(p = c(0), zero = c(0),n=c(0))

  p = ggplot(df, aes(x = p, y = zero,color=n))
  for( i in 0:6 ) {
    n = 2**i
    #  for( n in seq(1,8,by=1) ) {
    df = data.frame(p = Ps, zero = sapply(Ps, function(prob) {X = 10**6*rbinom(10000,n,prob)/rbinom(1000,N,prob); sd(X)/mean(X)}),
                    n = rep(format(n,width=3),length(Ps)))
    p = p + geom_line(data=df,size=2)
  }
  p = p + scale_color_brewer(palette = "Oranges")
  p = p+theme_minimal()
  p = p + labs(title="TPM coeff. of variation given capture", x="Capture Efficiency", y = "CV(TPM)")+scale_y_continuous(limits=c(0,2))
  
  
}