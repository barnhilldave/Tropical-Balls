# Sampling and Volume estimation
set.seed(123)
P<-matrix(c(0,-2,3,0,-2,5,0,1,0,0,2,2),4,3,TRUE)

B_p<-min_enc_ball(P)
B_r<-trop_bal.vert(B_p[[1]],B_p[[2]],3)
N<-1000
har_points<-matrix(0,N,ncol(P))
har_points1<-matrix(NA,0,ncol(P))
count<-0
x0<-c(0,0,0)
for (i in 1:N){
  x<-TropicalPolytope.extrapolation.HAR_v4(B_r,x0,I=50,k=2)
  har_points[i,]<-x
  px<-project_pi(P,x)
  d_p<-trop.dist(px,x)
  if (d_p<=1e-8){
    har_points1<-rbind(har_points1,x)
    count<-count+1
  }
  x0<-x
}
plot(har_points[,2],har_points[,3],pch=19,cex=.5,asp=1,col='darkgrey')
points(har_points1[,2],har_points1[,3],pch=19,cex=.5)

N<-10000
har_points<-matrix(0,N,ncol(P))
har_points1<-matrix(NA,0,ncol(P))
count1<-0
x0<-c(0,0,0)
for (i in 1:N){
  x<-TropicalPolytope.extrapolation.HAR_v4(B_r,x0,I=50,k=2)
  har_points[i,]<-x
  px<-project_pi(P,x)
  d_p<-trop.dist(px,x)
  if (d_p<=1e-8){
    har_points1<-rbind(har_points1,x)
    count1<-count+1
  }
  x0<-x
}
plot(har_points[,2],har_points[,3],pch=19,cex=.5,asp=1,col='darkgrey')
points(har_points1[,2],har_points1[,3],pch=19,cex=.5)

## Volume Estimation for N=1000
e<-dim(B_r)[2]
Vol_B_r<-e*B_p[[2]]^(e-1)
Vol_P<-count/1000*Vol_B_r

## Volume Estimation for N=1000
e<-dim(B_r)[2]
Vol_B_r<-e*B_p[[2]]^(e-1)
Vol_P<-count1/10000*Vol_B_r

# Volume using R functions
x0<-c(0,0,0)
Trop_Volume(B_r,P,x0,1000,50,d=ncol(P),B_p[[2]])

