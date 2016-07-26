#temp file

# A is a great Ab
sA = 0.1
bA = 0.001

# B is lousy Ab
sB = 0.05
bB = 0.005

A = matrix(c(sA, sA, bA, bA, 
             sB, bB, sB, bB, 
             sA*sA, sA*sA, bA*bA, bA*bA, 
             sB*sB, bB*bB, sB*sB, bB*bB, 
             sA*sB, sA*bB, sA*bB, bA*bB), ncol = 4, nrow = 5, byrow = "TRUE")





M = cbind( c( 1, 0, 0, 0), c(0,0,0,1), c(0.5, 0, 0, .5), c(0,0.5,0.5,0),c(.1,0,0,.9),c(.1,.2,.2,.5))

Y = A %*% M

