source("THAR_Funcs.R")

####### Rounding example to fit a minimum enclosing ball
####### around only the 3-trunk of a polytope that is not a polytrope

P<-matrix(c(0,0,0,0,0,1,3,1,0,1,2,5,0,2,5,10),4,4,TRUE)

# View the original polytope
draw.tpolytope.3d(P,'darkgrey')

B_r<-min_enc_ball(P) # minimum enclosing ball
B_ro<-trop_bal.vert(B_r[[1]],B_r[[2]],4) # Generating set for the ball

B_R<-max_ins_ball(P) # maximum inscribed ball

# Consider acceptance rate lower bound of the acceptance rate

ar<-(B_R[[1]]/B_r[[2]])^(4-1)
ar # .001

# Goal is to increase the acceptance rate by removing all e-2 tentacles
# The function converts the vertex representation of the tropical polytope
# into a hyperplane representation.  It then applies a double description
# algorithm to identify the pseudo-vertices of the 3-trunk.

PP<-rounding(P) # set of pseudo vertices

B_r_dd<-min_enc_ball(PP)
B_rm<-trop_bal.vert(B_r_dd[[1]],B_r_dd[[2]],4) # Generating set for the ball
# Size of the maximum inscribed ball does not change.
# Acceptance rate changes by over an order of magnitude.

ar_d<-(B_R[[1]]/B_r_dd[[2]])^(4-1) # 0.015625

# Compare a Volume estimate using the original minimum enclosing ball
# and the new.

# Before rou
N<-2000
x0<-c(0,0,0,0)
har_points<-matrix(0,N,4,TRUE)
har_points1<-matrix(NA,0,4,TRUE)
count=0
# Before rounding
for (i in 1:N){
  x<-TropicalPolytope.extrapolation.HAR_v4(B_ro,x0,I=50,k=3)
  har_points[i,]<-x
  px<-project_pi(P,x)
  d_p<-trop.dist(px,x)
  if (d_p<=1e-8){
    har_points1<-rbind(har_points1,x)
    count<-count+1
  }
  x0<-x
}

Vol_Br<-ncol(B_ro)*B_r[[2]]^(ncol(B_ro)-1)
p<-count/N
Vol_P<-p*Vol_Br

# After rounding
x0<-c(0,0,0,0)
har_points2<-matrix(0,N,4,TRUE)
har_points3<-matrix(NA,0,4,TRUE)
count1<-0
for (i in 1:N){
  x<-TropicalPolytope.extrapolation.HAR_v4(B_rm,x0,I=50,k=3)
  har_points2[i,]<-x
  px<-project_pi(P,x)
  d_p<-trop.dist(px,x)
  if (d_p<=1e-8){
    har_points3<-rbind(har_points1,x)
    count1<-count1+1
  }
  x0<-x
}

Vol_Brm<-ncol(B_rm)*B_r_dd[[2]]^(ncol(B_rm)-1)
p1<-count1/N
Vol_P_dd<-p1*Vol_Brm

Vol_P_dd # Volume estimate after rounding

