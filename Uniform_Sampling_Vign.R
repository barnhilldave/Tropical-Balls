##################################################################
# This code illustrates how to sample from a tropical polytope.  #
# that is not classically convex and also not a tropical simplex.#
##################################################################
source("THAR_funcs.R")

set.seed(123)
P<-matrix(c(0,-2,3,0,-2,5,0,1,0,0,2,2),4,3,TRUE)

A_P<-t(combn(seq(1,nrow(P)),dim(P)[2]))

B_r<-min_enc_ball(P)
d<-dim(P)[2]
B_p<-trop_bal.vert(B_r[[1]],B_r[[2]],3)
counts<-c(rep(0,nrow(A_P)))
C<-0
N<-2000
har_points<-matrix(0,N,dim(P)[2],TRUE)
har_points1<-matrix(NA,0,dim(P)[2],TRUE)
x0<-c(0,0,0)
for (i in 1:N){
  x<-TropicalPolytope.extrapolation.HAR_v4(B_p,x0,I = 50,k=2)
  har_points[i,]<-x
  x0<-x
  projp<-project_pi(P,x)
  dP<-trop.dist(projp,x)
  if(dP<=1e-8){
    C<-C+1
    har_points1<-rbind(har_points1,x)
  }
  for (j in 1:nrow(A_P)){
    proj<-project_pi(P[A_P[j,],],x)
    dp<-trop.dist(proj,x)
    if(dp<=1e-8){
      counts[j]<-counts[j]+1
    }
  }
}
ps<-counts/C
k=2
while (k<=nrow(P)){
  combos<-combinations(nrow(P),k)
  sums<-rep(0,nrow(combos))
  for(j in 1:length(sums)){
    sums[j]<-sum(counts[combos[j,1]],counts[combos[j,2]])
  }
  if(any(sums==C)){
    ind<-which(sums==C)
    inds<-combos[ind,]
    break
  }
  else{
    k=k+1
  }
}

p1<-ps[inds[1,1]]
p2<-ps[inds[1,2]]
P1<-P[A_P[inds[1,1],],]
P2<-P[A_P[inds[1,2],],]
M<-2000
har_points_p<-matrix(0,M,dim(P)[2],TRUE)
for (i in 1:M){
  p<-runif(1,0,1)
  if(p<=p1){
    x<-TropicalPolytope.extrapolation.HAR_v4(P1,x0,50,dim(P)[2]-1)
  }
  else{
    x<-TropicalPolytope.extrapolation.HAR_v4(P2,x0,50,dim(P)[2]-1)
  }
  har_points_p[i,]<-x
  x0<-x
}
plot(har_points_p[,2],har_points_p[,3],pch=19,col='black',
     cex=.5,asp=1,xlab='x2',ylab='x3',xlim=c(-2,2))

############### 3D example
# Polytope in R^4/R1
P<-matrix(c(0,0,2,0,0,0,0,0,0,4,-10,0,0,3,-3,5,0,4,2,10),5,4,TRUE)

# Find combinations of vertices that define all the simplices
A_P<-t(combn(seq(1,nrow(P)),dim(P)[2]))

# Calculate the minimum enclosing ball and define its generating set of vertices
B_r<-min_enc_ball(P)
d<-dim(P)[2]
B_p<-trop_bal.vert(B_r[[1]],B_r[[2]],4)

# Sample from the minimum enclosing ball and identify when points
# fall in each of the tropical simplices
counts<-c(rep(0,nrow(A_P)))
C<-0
N<-3000
har_points<-matrix(0,N,dim(P)[2],TRUE)
har_points1<-matrix(NA,0,dim(P)[2],TRUE)
x0<-c(0,0,0,0)
for (i in 1:N){
  x<-TropicalPolytope.extrapolation.HAR_v4(B_p,x0,I = 50,k=3)
  har_points[i,]<-x
  x0<-x
  projp<-project_pi(P,x)
  dP<-trop.dist(projp,x)
  if(dP<=1e-8){
    C<-C+1
    har_points1<-rbind(har_points1,x)
  }
  for (j in 1:nrow(A_P)){
    proj<-project_pi(P[A_P[j,],],x)
    dp<-trop.dist(proj,x)
    if(dp<=1e-8){
      counts[j]<-counts[j]+1
    }
  }
}
ps<-counts/C

# Considering only those points that fall in P, find the collection
# of tropical simplices that define the 3-trunk based on the
# proportion of points that fall in each simplex.
k=2
while (k<=nrow(P)){
  combos<-combinations(nrow(P),k)
  sums<-rep(0,nrow(combos))
  for(j in 1:length(sums)){
    sums[j]<-sum(counts[combos[j,1]],counts[combos[j,2]])
  }
  if(any(sums==C)){
    ind<-which(sums==C)
    inds<-combos[ind,]
    break
  }
  else{
    k=k+1
  }
}

# Here we find two of the simplices that define the 3-trunk
# We sample each according to the proportion of points 
# that fell in each simplex.
p1<-ps[inds[1]]
p2<-ps[inds[2]]
P1<-P[A_P[inds[1],],]
P2<-P[A_P[inds[2],],]
M<-3000
har_points_p<-matrix(0,M,dim(P)[2],TRUE)
for (i in 1:M){
  p<-runif(1,0,1)
  if(p<=p1){
    x<-TropicalPolytope.extrapolation.HAR_v4(P1,x0,50,dim(P)[2]-1)
  }
  else{
    x<-TropicalPolytope.extrapolation.HAR_v4(P2,x0,50,dim(P)[2]-1)
  }
  har_points_p[i,]<-x
  x0<-x
}

### Plot
plot3d(har_points_p[,2],har_points_p[,3],har_points_p[,4],
       xlab='x2',ylab='x3',zlab='x4',pch=19,col='red')

# With tropical polytope overlay
draw.tpolytope.3d(P,'lightgrey')
points3d(har_points_p[,2],har_points_p[,3],har_points_p[,4],
         xlab='x2',ylab='x3',zlab='x4',pch=19,col='red')
