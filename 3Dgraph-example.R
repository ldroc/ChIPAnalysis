library(plot3D)
N = 0
fname <- function() { N <<- N+1; paste0("~/Temp/pic",N,".png")}


M <- mesh(seq(0,1,length.out = 25),seq(0,1,length.out = 25))
P = M$x
Q = M$y

Th = 30
Ph = 10
mod = ComputePairwiseExpectations(s,a1,a2)
N1 = ComputeMarginalDistribution(s,a1)
N2 = ComputeMarginalDistribution(s,a2)
I = mod$V.obs < 1 
scatter3D(N1[I],N2[I], mod$V.obs[I],
          cex=.2,d = 1000,theta = Th,phi=Ph, type="d", xlim=c(0,1), ylim=c(0,1),zlim=c(0,1) )

#surf3D(P,Q,P*Q, d=1000,
#       border = "gray", box = T, bty = "b", panel.first = T,facets=T,theta = Th, phi=Ph)

surf3D(P,Q,pmin(P,Q), d=1000,
       border = "black", box = T, bty = "b", panel.first = T,facets=F, theta = Th, phi=Ph, add = TRUE)

surf3D(P,Q,pmax(0*P,P+Q-1), d=1000,
       border = "black", box = T, bty = "b", panel.first = T,facets=F,theta = Th, phi=Ph, add = TRUE)





surf3D(P,Q,pmin(P,Q) - P*Q, d=1000,
       border = "gray", box = T, bty = "b", panel.first = T,facets=T,zlim=c(-.25,.25),theta = Th, phi=Ph/2)

surf3D(P,Q,pmax(P*0,P+Q-1) - P*Q, d=1000,
       border = "gray", add=T, box = T, bty = "b", panel.first = T,facets=T,theta = Th, phi=Ph)

if(0) {
u.min = -0.5*pi
u.max = 0.5*pi
v.min = -pi
v.max = pi
M <- mesh(seq(u.min,u.max, length.out = 80),
          +           seq(v.min,v.max, length.out = 40))
u <- M$x
v <- M$y
x = cos(u)*cos(v) + 3*cos(u)*(1.5+sin(1.5*u)/2)
y = sin(u)*cos(v) + 3*sin(u)*(1.5+sin(1.5*u)/2)
z = sin(v)+2*cos(1.5*u)

png(fname())
surf3D(x, y, z, colvar = NA, col = "gray25", colkey = FALSE, xlim=c(-3,8), ylim=c(-7,9), zlim=c(-4,4),d=1000,
       border = "gray", box = T, bty = "b", panel.first = T,facets=F,theta = -40)
dev.off()
png(fname())
surf3D(x, y, z, colvar = NA, col = "gray25", colkey = FALSE, xlim=c(-3,8), ylim=c(-7,9), zlim=c(-4,4),d=1000,
       border = "gray", box = T, bty = "b", panel.first = T,facets=F,theta = 30, phi=20)
dev.off()
if(0){
  png(fname())
  surf3D(x, y, z, colvar = NA, col = "gray", colkey = FALSE, xlim=c(-3,8), ylim=c(-7,9), zlim=c(-4,4),
       border = "gray", facet=T, box = T, bty = "b", panel.first = T,d=1000,
       theta = 90,phi=0)
  dev.off()
  png(fname())
  surf3D(x, y, z, colvar = NA, col = "gray", colkey = FALSE, xlim=c(-3,8), ylim=c(-7,9), zlim=c(-4,4),
       border = "gray", facet=T, box = T, bty = "b", panel.first = T,d = 1000,
       theta = 180,phi=-90)
  dev.off()
  png(fname())
  surf3D(x, y, z, colvar = NA, col = "gray", colkey = FALSE, xlim=c(-3,8), ylim=c(-7,9), zlim=c(-4,4),
       border = "gray", facet=T, box = T, bty = "b", panel.first = T,d = 1000,
       theta = 0,phi=0)
  dev.off()
}
sample.u <- function(n) {
  u = rnorm(10*n)
  u = u[ u > u.min & u < u.max]
  u = u[1:n]
}
sample.v <- function(n) {
  v = rnorm(10*n)
  v = v[ v > v.min & v < v.max]
  v = v[1:n]
}


u = sample.u(100)
v = sample.v(100)
x = cos(u)*cos(v) + 3*cos(u)*(1.5+sin(1.5*u)/2)
y = sin(u)*cos(v) + 3*sin(u)*(1.5+sin(1.5*u)/2)
z = sin(v)+2*cos(1.5*u)

png(fname())
scatter3D(x,y,z,type="d",colkey=FALSE,xlim=c(-5,7), ylim=c(-7,9), zlim=c(-4,4),d=1000,
          cex=.2,theta = -40)
dev.off()
png(fname())
scatter3D(x,y,z,type="d",colkey=FALSE,xlim=c(-5,7), ylim=c(-7,9), zlim=c(-4,4),d=1000,
          cex=.2,theta = 30, phi = 20)
dev.off()


u = sample.u(500)
v = sample.v(500)
x = cos(u)*cos(v) + 3*cos(u)*(1.5+sin(1.5*u)/2)
y = sin(u)*cos(v) + 3*sin(u)*(1.5+sin(1.5*u)/2)
z = sin(v)+2*cos(1.5*u)

png(fname())
scatter3D(x,y,z,type="d",colkey=FALSE, xlim=c(-3,8), ylim=c(-7,9), zlim=c(-4,4),cex=.2, d=1000,theta = -40)
dev.off()
png(fname())
scatter3D(x,y,z,type="d",colkey=FALSE, xlim=c(-3,8), ylim=c(-7,9), zlim=c(-4,4),cex=.2, d=1000, phi = 20, theta = 30)
dev.off()

u = sample.u(10000)
v = sample.v(10000)
x = cos(u)*cos(v) + 3*cos(u)*(1.5+sin(1.5*u)/2)
y = sin(u)*cos(v) + 3*sin(u)*(1.5+sin(1.5*u)/2)
z = sin(v)+2*cos(1.5*u)

png(fname())
scatter3D(x,y,z,type="d",colkey=FALSE, xlim=c(-3,8), ylim=c(-7,9), zlim=c(-4,4),cex=.1, d = 1000,theta = 90,phi=0)
dev.off()
png(fname())
scatter3D(x,y,z,type="d",colkey=FALSE, xlim=c(-3,8), ylim=c(-7,9), zlim=c(-4,4),cex=.1, d = 1000,theta = 180,phi=-90)
dev.off()
png(fname())
scatter3D(x,y,z,type="d",colkey=FALSE, xlim=c(-3,8), ylim=c(-7,9), zlim=c(-4,4),cex=.1, d = 1000,theta = 0,phi=0)
dev.off()

u = sample.u(100000)
v = sample.v(100000)
x = cos(u)*cos(v) + 3*cos(u)*(1.5+sin(1.5*u)/2)
y = sin(u)*cos(v) + 3*sin(u)*(1.5+sin(1.5*u)/2)
z = sin(v)+2*cos(1.5*u)

png(fname())
hist(x,breaks="FD",freq=F,main="",xlab="",ylab="")
dev.off()
png(fname())
hist(y,breaks="FD",freq=F,main="",xlab="",ylab="")
dev.off()
png(fname())
hist(z,breaks="FD",freq=F,main="",xlab="",ylab="")
dev.off()

}