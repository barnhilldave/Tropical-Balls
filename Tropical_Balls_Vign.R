library(combinat)

# Construct maximum inscribed tropical ball for a polytrope in 2D
P<-matrix(c(0,0,0,0,2,5,0,3,1),3,3, TRUE)
v<-max_ins_ball(P)
plot(P[,2],P[,3],col='white',asp=1,xlab = 'x2',ylab='x3')

segments(0,0,0,3)
segments(2,5,0,3)
segments(2,5,3,5)
segments(3,5,3,1)
segments(3,1,2,0)
segments(0,0,2,0)
points(P[,2],P[,3],col='darkgrey',pch=19)

Trop_ball(v[[2]],v[[1]],cls='darkgrey',plt=FALSE)

# Construct minimum enclosing tropical ball for a polytrope in 2D

vv<-min_enc_ball(P)
Trop_ball(vv[[1]],vv[[2]],cls='white',plt=TRUE)
plot(P[,2],P[,3],col='white',asp=1,xlab = 'x2',ylab='x3')

segments(0,0,0,3)
segments(2,5,0,3)
segments(2,5,3,5)
segments(3,5,3,1)
segments(3,1,2,0)
segments(0,0,2,0)
points(P[,2],P[,3],col='darkgrey',pch=19)

# Construct maximum inscribed tropical ball for a polytrope in 2D

P<-matrix(c(0,0,0,0,2,5,0,3,1),3,3, TRUE)
v<-max_ins_ball(P)
plot(P[,2],P[,3],col='white',asp=1,xlab = 'x2',ylab='x3')

segments(0,0,0,3)
segments(2,5,0,3)
segments(2,5,3,5)
segments(3,5,3,1)
segments(3,1,2,0)
segments(0,0,2,0)
points(P[,2],P[,3],col='darkgrey',pch=19)

Trop_ball(v[[2]],v[[1]],cls='darkgrey',plt=FALSE)

# Construct maximum inscribed tropical ball for a generic tropical polytope in 2D
P<-matrix(c(0,-2,3,0,-2,5,0,1,0,0,2,2),4,3,TRUE)

# Identify all tropical simplices
A_P<-t(combn(seq(1,4),dim(P)[2])) # All combinations of row indices for the vertex set
R<-c()
for (i in 1:nrow(A_P)){
  K<-P[A_P[i,],]
  v<-max_ins_ball(K)
  R<-append(R,v[[1]])
}

R_ind<-which.max(R)
P_prime<-P[A_P[R_ind,],]

vv<-max_ins_ball(P_prime)
plot(P[,2],P[,3],col='white',asp=1,xlab = 'x2',ylab='x3')
segments(-2,3,-2,5)
segments(-2,5,2,5)
segments(2,5,2,2)
segments(2,2,1,1)
segments(1,0,1,3)
segments(-2,3,1,3)

Trop_ball(vv[[2]],vv[[1]],cls='gray',plt=FALSE)

points(P[,2],P[,3],col='darkgrey',pch=19)

# Construct minimum encompassing tropical ball for a generic tropical polytope that is not a
# tropical simplex in 2D
P<-matrix(c(0,-2,3,0,-2,5,0,1,0,0,2,2),4,3,TRUE)

vv<-min_enc_ball(P)
Trop_ball(vv[[1]],vv[[2]],cls='white',plt=TRUE)
segments(-2,3,-2,5)
segments(-2,5,2,5)
segments(2,5,2,2)
segments(2,2,1,1)
segments(1,0,1,3)
segments(-2,3,1,3)

points(P[,2],P[,3],col='darkgrey',pch=19)


# Construct maximum inscribed tropical ball for a generic tropical simplex in 3D
P<-matrix(c(0,0,0,0,0,1,3,1,0,1,2,5,0,2,5,10),4,4,TRUE)

vv<-max_ins_ball(P)
draw.tpolytope.3d(P,'lightgrey')
Trop_ball(vv[[2]],vv[[1]],cls='black',plt=FALSE)


