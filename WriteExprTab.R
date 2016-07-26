library(ctc)

WriteExpressionTab <- function( X, file, Genome ) {
  n = ncol(X)+2
  Y = matrix(ncol = n, nrow = nrow(X))
  Y[,3:n] = X
  colnames(Y) = c("","", colnames(X))
  Y[,1] = rownames(X)
  Y[,2] = AccToGeneName(rownames(X), Genome)
  
  r2cluster(Y, file=file, description=TRUE, labels=TRUE)
}

