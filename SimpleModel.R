

## parameters for ab 
##
# X is in [AB, Ab, aB, ab]
#
# Parameters  sA, bA, sB, bB, Ma, Mb
#
# output  A, B, A:B, B:A, A:A, B:B
#

abModelBasicCounts <- function( X, p ) {
  # ugly code
  Ma = p[["Ma"]]
  Mb = p[["Mb"]]
  sA = p[["sA"]]
  sB = p[["sB"]]
  bA = p[["bA"]]
  bB = p[["bB"]]
  XAB = X[["AB"]]
  XAb = X[["Ab"]]
  XaB = X[["aB"]]
  Xab = X[["ab"]]
  
  unlist(list ( 
    A =  Ma * ( sA * (XAB+XAb) + bA*(XaB+Xab) ),
    B =  Mb * ( sB * (XAB+XaB) + bB*(XAb+Xab) ),
    AB = Ma * ( sA*sB*XAB + sA*bB*XAb + bA*sB*XaB + bA*bB*Xab ),
    BA = Mb * ( sA*sB*XAB + sA*bB*XAb + bA*sB*XaB + bA*bB*Xab ), 
    AA = Ma * ( sA * sA * (XAB+XAb) + bA * bA * (XaB+Xab) ),
    BB = Mb * ( sB * sB * (XAB+XaB) + bB * bB * (XAb+Xab) )
  ))
}

abModelBasicCountsGrad <- function( X, p ) {
  # ugly code
  Ma = p[["Ma"]]
  Mb = p[["Mb"]]
  sA = p[["sA"]]
  sB = p[["sB"]]
  bA = p[["bA"]]
  bB = p[["bB"]]
  XAB = X[["AB"]]
  XAb = X[["Ab"]]
  XaB = X[["aB"]]
  Xab = X[["ab"]]
  
  
  matrix( c(Ma * sA, Ma * sA, Ma * bA, Ma * bA, sA * (XAB+XAb) + bA*(XaB+Xab), 0, Ma * (XAB+XAb), 0, Ma*(XaB+Xab), 0 ,  # A
            Mb * sB, Mb * bB, Mb * sB, Mb * bB, 0, ( sB * (XAB+XaB) + bB*(XAb+Xab) ), 0, Mb*(XAB+XaB), 0, Mb*(XAb+Xab), # B
            Ma*sA*sB, Ma*sA*bB, Ma*bA*sB, Ma*bA*bB, sA*sB*XAB + sA*bB*XAb + bA*sB*XaB + bA*bB*Xab, 0, #AB
            Ma*(sB*XAB + bB*XAb), Ma*(sA*XAB + bA*XaB), Ma*(sB*XaB + bB*Xab), Ma*(sA*XAb + bA*Xab),
            Mb*sA*sB, Mb*sA*bB, Mb*bA*sB, Mb*bA*bB, 0, sA*sB*XAB + sA*bB*XAb + bA*sB*XaB + bA*bB*Xab, #BA
            Mb*(sB*XAB + bB*XAb), Mb*(sA*XAB + bA*XaB), Mb*(sB*XaB + bB*Xab), Mb*(sA*XAb + bA*Xab),
            Ma*sA*sA, Ma*sA*sA, Ma*bA*bA, Ma*bA*bA,sA*sA*(XAB+XAb) + bA*bA*(XaB+Xab), 0, 2*Ma*sA*(XAB+XAb), 0, 2*Ma*bA*(XaB+Xab), 0, #AA
            Mb*sB*sB, Mb*bB*bB, Mb*sB*sB, Mb*bB*bB, 0, sB * sB*(XAB+XaB) + bB * bB*(XAb+Xab), 0, 2*Mb*sB*(XAB+XaB),0, 2*Mb*bB*(XAb+Xab) 
          ), nr = 6, nc = 10, byrow = T,   
          dimnames = list( c("A", "B", "AB", "BA", "AA", "BB"),
                           c("AB", "Ab", "aB", "ab", "Ma", "Mb", "sA", "sB", "bA", "bB"))
  )
}

abModelQuadLoss <- function( N, Y ) {
  return(-sum( (N-Y)**2 ))
}

abModelPoissonLoss <- function( N, Y ) {
  return(-sum( N*log(Y) - Y ) )
}

abModelPoissonLossGrad <- function( N, Y, G ) {
  -colSums((N/Y-1)*G)
}

abModelKOCounts <- function( X, p, Mda, Mdb ) {
  Y = abModelBasicCounts( X, p )
  Xda = X
  Xda["AB"] = 0
  Xda["Ab"] = 0
  pda = p
  pda["Ma"] = Mda[["Ma"]]
  pda["Mb"] = Mda[["Mb"]]
  Yda = abModelBasicCounts( Xda, pda )
  Xdb = X
  Xdb["AB"] = 0
  Xdb["aB"] = 0
  pdb = p
  pdb["Ma"] = Mdb[["Ma"]]
  pdb["Mb"] = Mdb[["Mb"]]
  Ydb = abModelBasicCounts( Xdb, pdb )
  c(Y,Yda,Ydb)
}

abModelKOCountsGrad <- function( X, p, Mda, Mdb ) {
  G = abModelBasicCountsGrad( X, p )
  Z = matrix(0, nr = 6, nc = 2)
  G = cbind( G, Z, Z)

    Xda = X
  Xda["AB"] = 0
  Xda["Ab"] = 0
  pda = p
  pda["Ma"] = Mda[["Ma"]]
  pda["Mb"] = Mda[["Mb"]]
  Gda = abModelBasicCountsGrad( Xda, pda )
  Gda = cbind( Gda[,1:4], Z, Gda[,7:10], Gda[,5:6], Z)
  Gda[,"AB"] = 0
  Gda[,"Ab"] = 0
  Xdb = X
  Xdb["AB"] = 0
  Xdb["aB"] = 0
  pdb = p
  pdb["Ma"] = Mdb[["Ma"]]
  pdb["Mb"] = Mdb[["Mb"]]
  Gdb = abModelBasicCountsGrad( Xdb, pdb )
  Gdb = cbind( Gdb[,1:4], Z, Gdb[,7:10], Z, Gdb[,5:6])
  Gdb[,"AB"] = 0
  Gdb[,"aB"] = 0
  rbind(G,Gda,Gdb)
}

MaskAbPref = c(rep(1,6),rep(0,4), rep(1,4))
MaskX = c(rep(0,4),rep(1,10))

abModelKOOptimize <- function ( N, 
                                X = list( AB = .25, Ab = .25, aB = .25, ab = .25), 
                                p = list( Ma = 1000000, Mb = 1000000, sA = .1, sB = .1, bA = 0.001, bB = 0.001), 
                                Mda = list( Ma = 1000000, Mb = 1000000), 
                                Mdb = list( Ma = 1000000, Mb = 1000000 ),
                                Mask = rep(1,14)
) {
  x = c(X, p, Mda, Mdb )
  
  loss <- function( x ) {
    X[1:4] = x[1:4]
    p[1:6] = x[5:10]
    Mda[1:2] = x[11:12]
    Mdb[1:2] = x[13:14]
    Y = abModelKOCounts( X, p, Mda, Mdb )
    abModelPoissonLoss(N,Y)
  }
  
  grad <- function( x ) {
    X[1:4] = x[1:4]
    p[1:6] = x[5:10]
    Mda[1:2] = x[11:12]
    Mdb[1:2] = x[13:14]
    Y = abModelKOCounts( X, p, Mda, Mdb )
    G = abModelKOCountsGrad( X, p, Mda, Mdb )
    
    abModelPoissonLossGrad(N,Y,G)*Mask
  }

  up = c( rep(1,4), rep(Inf,2), rep(1,4), rep(Inf,4))
  opt = optim(x, fn = loss, gr = grad, lower = 1e-5, upper = up, method = "L-BFGS-B")
  print(opt)
  x = opt$par
  X[1:4] = x[1:4]
  p[1:6] = x[5:10]
  Mda[1:2] = x[11:12]
  Mdb[1:2] = x[13:14]
  return(list(X = X, p = p, Mda = Mda, Mdb = Mdb))
}


if( 0 ) {
  # test code
  p = list( Ma = 1000, Mb = 500, sA = 0.1, sB = 0.05, bA = 0.001, bB = 0.005 )
  X = list( AB = .2, Ab = .1, aB = .3, ab = .2)
  Mda = list( Ma = 1000, Mb = 1000)
  Mdb = list( Ma = 1000, Mb = 1000)
  
  X = list( AB = 1, Ab = 0, aB = 0, ab = 0)
  
  cnts = read.csv("~/Data/CoChIP/HisMut-K4/HisMut/HisMut-all-counts.csv", row.names = 1)

  BuildStrainNameList <- function ( s, A, B ) {
    a = paste0("H3",A)
    b = paste0("H3", B)
    paste0(s,"-", c( paste0(a,"-Input"), paste0(b,"-Input"), 
                     paste0(a,"-", b), paste0(b,"-", a),
                     paste0(a,"-", a), paste0(b,"-", b) ) )
  }
  
  BuildNameList <- function( s, A, B ) {
    do.call("c", lapply(s, function(x) BuildStrainNameList( x, A, B) ))
  }
  
  N.K4.K18 =  cnts[BuildNameList(c("WT", "K4R", "K18R"), "K4me3", "K18ac"),]
  N.K4.K36 =  cnts[BuildNameList(c("WT", "K4R", "K36R"), "K4me3", "K36me3"),]
}

testGrad <- function( X, p ) {
  eps = 1e-7
  Y = abModelBasicCounts( X, p )
  G = abModelBasicCountsGrad( X, p )
  Xt = X
  pt = p
  g = abModelPoissonLossGrad( N.K4.K18[1:6], Y, G)
  l = abModelPoissonLoss(N.K4.K18[1:6], Y)
  Gt = NULL
  gt = NULL
  for( i in 1:length(X) ) {
    Xt = X
    Xt[i] = Xt[[i]] + eps
    Yt = abModelBasicCounts( Xt, pt)
    lt = abModelPoissonLoss(N.K4.K18[1:6], Yt)
    
    if( is.null(Gt)) {
      Gt = (Yt - Y)/eps
      gt = c((lt-l)/eps)
    } else {
      Gt = cbind( Gt, (Yt - Y)/eps)
      gt = c( gt,(lt-l)/eps )
    }
  }
  Xt = X  
  for( i in 1:length(p) )
  {
    pt = p
    pt[i] = pt[[i]] + eps
    Yt = abModelBasicCounts( Xt, pt)
    lt = abModelPoissonLoss(N.K4.K18[1:6], Yt)
    Gt = cbind( Gt, (Yt - Y)/eps)
    gt = c( gt,(lt-l)/eps )
  }
  colnames(Gt) = colnames(G)
  names(gt) = names(g)
  print(Gt)
  print(G)
  print( G - Gt )
  print(gt)
  print(g)
  print(g - gt)
#  print(N.K4.K18[1:6]/Y-1)
#  print((N.K4.K18[1:6]/Y-1 )*(G-Gt))
  print(colSums((N.K4.K18[1:6]/Y-1 )*(G-Gt)))
  
}