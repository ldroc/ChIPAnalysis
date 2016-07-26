
ccSimulateSeq <- function(F, N) {
  rpois(length(F), N*F)
}

ccSimulateChIP <- function( par, Design, X ) {
  as.vector(t(par) %*% Design %*% X)
}

#
# T - HMM transition matrix
ccGenerateDetXs <- function( T, n ) {
  k = dim(T)[[1]]
  X = matrix( 0, ncol = n, nrow = k)
 

  X[,1] = rmultinom(1,1,rep(1,k))
  for(i in 2:n) {
    X[,i] = rmultinom(1,1,T %*% X[,i-1])
  }
  X
}

ccGenerateXs <- function( T, n, conv = c(0.01, 0.05,0.19,0.5,0.19,0.05, 0.01)) {
  m = length(conv)
  k = dim(T)[[1]]
  DX = ccGenerateDetXs(T,n+m)
  X = matrix( 0, ncol = n, nrow = k)
  for( i in 1:m )
    X = X+ (conv[i] * DX[,i:(n+i-1)])
  
  X
}

ccSimulateCounts = function( T, n, par, multiDesign, Ns) {
  k = length(multiDesign)
  Xs = ccGenerateXs( T, n)
  Counts = do.call(rbind,lapply(1:k, function(i) ccSimulateSeq( ccSimulateChIP(par, multiDesign[[i]], Xs), Ns[[i]])))
  list( Counts = Counts, Xs = t(Xs))
}

simulationT = matrix( c( 0.5, 0.0, 0.5, 0.0,  #both
                         0.3, 0.7, 0.0, 0.0,  #K4
                         0.0, 0.0, 0.8, 0.2,  #K36
                         0.0, 0.7, 0.0, 0.3  # none
                        ), ncol = 4, nrow = 4, byrow = TRUE )

