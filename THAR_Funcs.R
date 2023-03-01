#########################################################################
## Author :  Ruriko Yoshida and David Barnhill                         ##
## Date   :  March 26th 2022                                           ##
## Update :  Mark 1st 2023.                                            ##
## Program:  These functions are used for investigating                ##
##           characteristics of tropical polytopes and                 ##
##           data derivatives.                                         ##
#########################################################################
library(MASS)
library(ape)
library(phangorn)
library(lpSolve)
library(lpSolveAPI)
library("miscTools")
library(Matrix)
library(rcdd)
library( rgl )
library(magick)
library(gtools)
library(combinat)

printf <- function(...) invisible(print(sprintf(...)))



# This function is used for calculating R-Square
fermatweberdistance <- function(datamatrix) {
  n = dim(datamatrix)[1]
  print(n)
  m = dim(datamatrix)[2]
  lprec <- make.lp(0, n+m)
  objective = mat.or.vec(n+m,1)
  for (i in seq(n)) {
    objective[i] = 1
  }
  set.objfn(lprec, objective)
  for (i in seq(n)) {
    for (j in seq(m)) {
      for (k in seq(m)) {
        v = mat.or.vec(n+m,1)
        v[i] = 1
        v[n+k] = 1
        v[n+j] = -1
        add.constraint(lprec, v, ">=", datamatrix[i,k] - datamatrix[i,j])
      }
    }
  }
  solve(lprec)
  return((get.variables(lprec)))
}

trop.dist <- function(D1, D2){
  return(max(D2 - D1) - min(D2 - D1))
}

normaliz.tree <- function(D, h = 1){
  d <- length(D)
  a <- max(D) - h
  x <- D - rep(a, d)
  for(i in 1:d)
    if(x[i] < 0)
      x[i] <- 0
  return(x)
}

### The first coordinate must be 0.
normaliz.vector <- function(D){
  return(D - rep(D[1], length(D)))
}

normaliz.vectors <- function(DD){ ## DD is a matrix
  d <- dim(DD)
  D1 <- DD
  for(i in 1:d[1])
    D1[i, ] <- D1[i, ] - rep(D1[i, 1], d[2])
  return(D1)
}

normaliz.polytope <- function(M){
  d <- dim(M)
  for(i in 1:d[1])
    M[i, ] <- normaliz.vector(M[i,])
  return(M)
}


TLineSeg <- function(D1, D2){
  d <- length(D1)
  if(length(D1) != length(D2))
    warning("dimension is wrong!")
  index <- order(D2 - D1)
  lambda <- (D2 - D1)[index]
  x1 <- rep(0, d)
  segment <- list()
  for(j in 1:d){
    for(i in 1:d){
      x1[i] <- 0
      x1[i] <- max(lambda[j] + D1[i], D2[i])
    }
    segment[[j]] <- x1
  }
  return(segment)
}

TLineSeg_min <- function(D1, D2){
  d <- length(D1)
  if(length(D1) != length(D2))
    warning("dimension is wrong!")
  index <- order(D2 - D1,decreasing = TRUE)
  lambda <- (D2 - D1)[index]
  x1 <- rep(0, d)
  segment <- list()
  for(j in 1:d){
    for(i in 1:d){
      x1[i] <- 0
      x1[i] <- min(lambda[j] + D1[i], D2[i])
    }
    segment[[j]] <- normaliz.vector(x1)
  }
  return(segment)
}


Points.TLineSeg <- function(D1, D2, k = 20){
  d <- length(D1)
  x <- matrix(rep(0, d*k), k, d)
  if(length(D1) != length(D2))
    warning("dimension is wrong!")
  index <- order(D2 - D1)
  lambda <- (D2 - D1)[index]
  L <- lambda[d] - lambda[1]
  l <- runif(1, min=lambda[1], max=lambda[d])
  for(j in 1:(k-1))
    for(i in 1:d){
      x[j, i] <- max((lambda[1] + (j * L)/k) + D1[i], D2[i])
    }
  
  return(x)
}

HAR.TLineSeg <- function(D1, D2){
  d <- length(D1)
  if(length(D1) != length(D2))
    warning("dimension is wrong!")
  index <- order(D2 - D1)
  lambda <- (D2 - D1)[index]
  l <- runif(1, min=lambda[1], max=lambda[d])
  x1 <- rep(0, d)
  
  for(i in 1:d){
    x1[i] <- max(l + D1[i], D2[i])
  }
  
  return(x1)
}


#### HAR on a tropical polytope

project_pi<-function(D_s,D){
  d <- dim(D_s)
  lambda <- rep(0, d[1])
  for(i in 1:d[1])
    lambda[i] <- min(D - D_s[i,])
  
  x <- rep(0, d[2])
  for(i in 1:d[2])
    x[i] <- max((lambda+D_s)[, i])
  return(normaliz.vector(x))
}


sample.pi <- function(D_s, a, b){
  d <- dim(D_s)
  #index <- sample(1:d[1], 1)
  lambda <- rep(0, d[1])
  #lambda[index] <- runif(1, min=a, max=b)
  for(i in 1:d[1])
    lambda[i] <- runif(1, min=a[i], max=b[i])
  #for(i in 1:d[1])
  #    lambda[i] <- min(D - D_s[i,])
  
  x <- rep(0, d[2])
  for(i in 1:d[2])
    x[i] <- max((lambda+D_s)[, i])
  return(x)
}


### Input: rxd matrix and each row of the matrix is a vertex of the polytope.
Intervals <- function(D_s){
  d <- dim(D_s)
  L <- matrix(rep(0, d[1]*d[2]), d[2], d[1])
  for(i in 2:d[2])
    L[i, i] <- max(D_s[1,] - D_s[i,])
  for(j in 1:d[2]){
    for(k in 2:d[1]){
      if(k != j){
        L[j, k] <- L[j, j] - max(D_s[k,] - D_s[j,])
      }
    }
  }
  # print(L)
  bounds <- matrix(rep(0, 2*d[1]), 2, d[1])
  for(i in 1:d[1])
    bounds[1, i] <- min(L[,i])
  for(i in 1:d[1])
    bounds[2, i] <- max(L[,i])
  return(bounds)
}



#### HAR algorithm on the space of ultrametrics
### n is the number of leaves

Ultrametrics.HAR <- function(x0, n, I = 1, h = 1){
  d <- length(x0)
  x <- normaliz.tree(x0, h)
  a <- h
  
  for(i in 1:I){
    x1 <- normaliz.tree(x, h)
    D0 <- symMatrix(runif(choose((n+1), 2), 0, a), n)
    for(k in 1:n)
      D0[k, k] <- 0
    tree <- upgma(D0)
    #tree <- rcoal(n)
    tree$edge.length <- tree$edge.length/max(tree$edge.length)
    D <- cophenetic(tree)
    v <- D[lower.tri(t(D))] #rep(0, d)
    #t <- 1
    #for(ii in 1:(n-1))
    #    for(jj in (ii+1):n){
    #        v[t] <- D[ii, jj]
    #        t <- t + 1
    #    }
    x <- HAR.TLineSeg(normaliz.tree(x1, h), normaliz.tree(v, h))
  }
  
  return(normaliz.tree(x, h))
}




### D is the set of vertices for tropical polytope with row is a vertex
### Only for e = 4.
pre.draw.tpolytope.3d <- function(D, v,c){
  d <- dim(D)
  seg <- matrix(rep(0, 3*choose(d[1], 2)*(2*3)), 3*choose(d[1], 2), 6)
  counter <- 1
  for(i in 1:(d[1] - 1)){
    for(j in (i+1):d[1]){
      t <- TLineSeg(D[i, ], D[j, ])
      for(k in 1:3){
        seg[counter, 1:3] <- normaliz.vector(t[[k]])[2:4]
        seg[counter, 4:6] <- normaliz.vector(t[[k+1]])[2:4]
        counter <- counter + 1
      }
    }
  }
  
  segments3d(x=as.vector(t(seg[1:(counter-1), c(1,4)])),
             y=as.vector(t(seg[1:(counter-1), c(2,5)])),
             z=as.vector(t(seg[1:(counter-1), c(3,6)])), col = c, lwd = .2,tcl=-.9)
  for(i in 1:v)
    spheres3d(D[i, 2:4], r = 0.1, color = "black") 
  #rgl.points(D[i, 2:4], color = "blue", size = 10)
  
  axes3d()
  title3d(xlab="X",ylab="Y",zlab="Z")
}

draw.tpolytope.3d <- function(D,c){
  d <- dim(D)
  D1 <- D
  for(i in 1:(d[1] - 1)){
    for(j in (i+1):d[1]){
      M <- Points.TLineSeg(D[i, ], D[j, ])
      D1 <- rbind(D1, M)
    }
  }
  pre.draw.tpolytope.3d(D1, d[1],c)
  
}

###### Check whether x is in a tropical polytope P or not.

Check.onto.Tpoly <- function(D_s, x){
  ### D_s is a set of vertices for a tropical polytope.
  ### Each row of D_s is a vertex of the tropical polytope.
  ### We would like to check whether x is on the tropical polytope or not.
  ### If x is on the tropical polytope then it returns 1.  If not it returns 0.
  y <- project_pi(D_s, x)
  res <- 0
  if(trop.dist(x, y) == 0) res <- 1
  return(res)
}

##### Sampling from Gaussian
### mu: center of mass
### D: set of vertices of a tropical polytope
### s: standard deviation
### I: number of iterations
trop.Gaussian <- function(D, x0,mu, s, n){
  
  d <- dim(D)
  mu <- normaliz.vector(mu)
  N <- matrix(rep(0, n*d[2]), n, d[2])
  #x0 <- mu
  i <- 1
  while(i <= n){
    print(i)
    x1 <- TropicalPolytope.V.HAR(D, x0, I = 50, k = 3)
    #x1 <- TropicalPolytope.extrapolation.HAR_v4(D, x0, I = 50,k=3)
    x1 <- normaliz.vector(x1)
    (r <- exp(-trop.dist(mu, x1)^2/s)/exp(-trop.dist(mu, x0)^2/s))
    #r<-exp(-trop.dist(mu, x1)/s)/exp(-trop.dist(mu, x0)/s)
    if(runif(1) < r){
      x0 <- x1
      N[i, ] <- x0
      N[i, ] <-  normaliz.vector(N[i, ])
      i <- i + 1
    }
  }
  return(N)
}


###########################################
## Author :  David Barnhill
## Date   :  October 12th 2022
## Program:  This code produces a tropical ball in two or three dimensions
## Input  :  Central point, v; tropical radius of ball; color to use for lines or 
##           shading shading level (Defaults to black); of each side, a, (1 equals opaque); 
##           fill logical indicating whether to fill each facet with color. 
##           If fill is 'TRUE' then then the one dimensional facets (line segments) 
##           show (Default is FALSE); plot logical indicating whether to 
##           execute a new plot of the ball or to plot it on an existing 
##           plot (Default is TRUE); border for 2D plots (Default is black). 
##           For 2D plots, if a plot is not already open, plot must be TRUE.  
##           For 3D plots, new plot is not necessary.
##           Either way, if a plot is already present, set plot to 'FALSE' and
##           the ball will overlay on the current plot.
## Output :  Plot or overlay of Euclidean representaiton oftropical ball in 2D or 3D.
## Execute: type in R as 
##
#############################################

Trop_ball<-function(v,d,a=1,cls='black',cent.col='black',fil=TRUE,plt=TRUE,bord='black'){
  dm<-length(v)-1
  if(dm==2){
    if (plt==TRUE){
      plot(v[2],v[3],type='n',xlim =c(v[2]-d-.25,(v[2]+d+.25)) ,
           ylim =c(v[3]-d-.25,(v[3]+d+.25)) ,
           xlab='x2',ylab='x3')
    }
    polygon(c(v[2]-d,v[2]-d,v[2],v[2]+d,v[2]+d,v[2]),c(v[3]-d,v[3],v[3]+d,v[3]+d,v[3],v[3]-d),col=cls,density=70,border = bord,lwd=2)
    points(v[2],v[3],pch=19,col=cent.col)
  }
  
  if(dm==3){
    if (plt==TRUE){
      plot3d(v[2],v[3],v[4],pch=19,col=cent.col,size=6,xlim =c(v[2]-d-.25,(v[2]+d+.25)) ,
             ylim =c(v[3]-d-.25,(v[3]+d+.25)) ,
             zlim=c(v[4]-d-.25,(v[4]+d+.25)),
             xlab='x2',ylab='x3',zlab='x4')
    }
    
    polygon3d(c(v[2]-d,v[2]-d,v[2]-d,v[2]-d),c(v[3],v[3]-d,v[3]-d,v[3]),c(v[4],v[4],v[4]-d,v[4]-d),
              coords = c(2,3),fill = fil,alpha=a,col=cls) 
    polygon3d(c(v[2],v[2]-d,v[2]-d,v[2]),c(v[3]-d,v[3]-d,v[3]-d,v[3]-d),c(v[4],v[4],v[4]-d,v[4]-d),
              coords = c(1,3),fill = fil,alpha=a,col=cls)
    polygon3d(c(v[2],v[2]-d,v[2]-d,v[2]),c(v[3],v[3],v[3]-d,v[3]-d),c(v[4]-d,v[4]-d,v[4]-d,v[4]-d),
              coords = c(1,2),fill = fil,alpha=a,col=cls)
    polygon3d(c(v[2]+d,v[2]+d,v[2]+d,v[2]+d),c(v[3],v[3]+d,v[3]+d,v[3]),c(v[4],v[4],v[4]+d,v[4]+d),
              coords = c(2,3),fill = fil,alpha=a,col=cls)
    polygon3d(c(v[2],v[2]+d,v[2]+d,v[2]),c(v[3]+d,v[3]+d,v[3]+d,v[3]+d),c(v[4],v[4],v[4]+d,v[4]+d),
              coords = c(1,3),fill = fil,alpha=a,col=cls)
    polygon3d(c(v[2],v[2]+d,v[2]+d,v[2]),c(v[3],v[3],v[3]+d,v[3]+d),c(v[4]+d,v[4]+d,v[4]+d,v[4]+d),
              coords = c(1,2),fill = fil,alpha=a,col=cls)
    polygon3d(c(v[2]-d,v[2],v[2],v[2]-d),c(v[3],v[3]+d,v[3]+d,v[3]),c(v[4],v[4]+d,v[4],v[4]-d),
              coords = c(2,3),fill = fil,alpha=a,col=cls)
    polygon3d(c(v[2]-d,v[2],v[2],v[2]-d),c(v[3],v[3]+d,v[3],v[3]-d),c(v[4],v[4]+d,v[4]+d,v[4]),
              coords = c(2,3),fill = fil,alpha=a,col=cls)
    polygon3d(c(v[2]-d,v[2],v[2]+d,v[2]),c(v[3]-d,v[3],v[3],v[3]-d),c(v[4],v[4]+d,v[4]+d,v[4]),
              coords = c(1,2),fill = fil,alpha=a,col=cls)
    polygon3d(c(v[2],v[2]+d,v[2]+d,v[2]),c(v[3]-d,v[3],v[3],v[3]-d),c(v[4],v[4]+d,v[4],v[4]-d),
              coords = c(2,3),fill = fil,alpha=a,col=cls)
    polygon3d(c(v[2],v[2],v[2]+d,v[2]+d),c(v[3]-d,v[3],v[3]+d,v[3]),c(v[4]-d,v[4]-d,v[4],v[4]),
              coords = c(2,3),fill = fil,alpha=a,col=cls)
    polygon3d(c(v[2]-d,v[2],v[2]+d,v[2]),c(v[3],v[3],v[3]+d,v[3]+d),c(v[4]-d,v[4]-d,v[4],v[4]),
              coords = c(1,2),fill = fil,alpha=a,col=cls)
    if (plt==FALSE){
      points3d(v[2],v[3],v[4],pch=19,col='white',size=6)
    }
  }
}


#### Calculate distance from a point to a max hyperplane ####
trop.dist.hyp_max<-function(O,bp){
  x<-O+bp
  x_prime<-x[order(x,decreasing = TRUE)]
  dist<-x_prime[1]-x_prime[2]
  return(dist)
}

#### Calculate distance from a point to a min hyperplane ####
trop.dist.hyp_min<-function(O,bp){
  x<-O+bp
  x_prime<-x[order(x,decreasing = FALSE)]
  (dist<-x_prime[2]-x_prime[1])
  return(dist)
}

#### HAR Extrapolation with fixed nu. Used for polytropes ######
# Input: Vertices of polytope, D_s; starting point, x0; number #
# of steps, I; cardinality, k, of subset, U.                   #
# Output: Point x1 in polytope                                 #
################################################################

TropicalPolytope.extrapolation.HAR_v4 <- function(D_s, x0, I = 1,k=1){
  d <- dim(D_s) # dimension of matrix of vertices
  D <- normaliz.polytope(D_s) #D_s
  D_bp<-D 
  x <- normaliz.vector(x0) #normalize x0
  
  for(i in 1:I){
    x1 <- x
    v <- x
    u <- x
    ind<-seq(1,d[1]) 
    (index <- sample(ind, k)) # Randomly choose k vertices
    U <- D[index,] # Subset U
    V <- D[-index,] # Complement of U
    if(k == 1) # If subset k=1 then this is just a vertex
      u <- U   # the projection of x0 on a vertex is just the vertex
    else{
      u <- project_pi(U, x1) # If k>1 project x0 onto the polytope
    }                      # defined by the subset of vertices
    if(k == (d[1] - 1)) # If the cardinality of U is e-1 
      v <- V            # the complement is just the leftover vertex
    else
      v <- project_pi(V, x1) # Otherwise it is a projection
    Gam<-unique(TLineSeg(v,u)) # Get unique bendpoints in a tropical line
    bps<-normaliz.polytope(matrix(unlist(Gam),length(Gam),d[2],byrow=TRUE)) # Normalize the bendpoints
    # We calculate the tropical distance of the bendpoints (not the end points) of the tropical
    # line segment to each boundary vertex defining a hyperplane.
    # If there is more than one bendpoint, calculate the distance of each one.  If not then just
    # conduct HAR on the line segment.
    if(nrow(bps)>2){ 
      l<-matrix(u,1,ncol(bps),TRUE) # Matrix only consisting of the starting point.
      bp<-bps[2:(nrow(bps)),] # bendpoints of the line segment without the endpoints
      t=0
      while (t <nrow(bp)){
        t=t+1
        dists<-as.vector(apply(D,1,function(x) trop.dist.hyp_min(-x,bp[t,]))) # Calculates distance of bend point to all vertices
        if(all(dists>1e-8)){
          l<-rbind(l,bp[t,])
        }
        else{
          l<-rbind(l,bp[t,]) # If there is a zero add the point and then exit the loop.
          break
        }
      }
      x<-HAR.TLineSeg(l[1,], l[nrow(l),]) # Conduct HAR on line segment between u and the last bendpoint.
      x<-normaliz.vector(x)
    }
    else{
      x<-HAR.TLineSeg(u,v) # Conduct HAR if there are no bendpoints.
      x<-normaliz.vector(x)
    }
  }
  
  return(normaliz.vector(x))
}

####### Tropical Ball Vertices #################################################
# Given a point, x in R^e/R1 which representing the center of a tropical ball. #
# Calculate the vertices defining the tropical convex hull of the tropical.    # 
# ball.                                                                        #
# Input: x0 in R^e/R1                                                          #
# Output: R^(e x e) matrix representing tropical convex hull                   #
################################################################################
trop_bal.vert<-function(x,d,dm){
  A<-matrix(0,dm,dm,TRUE)
  i=1
  while (i < dm){
    i=i+1
    a<-c(rep(0,dm))
    a[i]<-d
    A[i,]<-x+a
  }
  A[1,]<-x-c(0,rep(d,(dm-1)))
  return(A)
}

####### Volume of a Tropical Polytope ##########################################
# Given a minimum enclosing tropical ball, B, with radius, R, surrounding a    #
# tropical polytope, P, sample from B and determine the proportion of points.  #
# that fall inside B and P.  Multiply this proportion by the volume of B to    #
# calculate the volume of P.
# Input: Minimum encompassing ball, B; convex hull of tropical polytope, P;    #
#        starting point x0; sample size, S; iterations for HAR algorithm,i;    #
#        dimension of polytope P in R^d/R1, d; radius of B, R.                 #
# Output (list): ratio of points in P to S; Volume of the tropical ball, B;    #
#                Volume of P; sampled points (matrix); sampled points in P     #
################################################################################
Trop_Volume<-function(B,P,x0,S,i,d=ncol(P),R){
  count<-0
  for (i in (1:S)){
    print(i)
    x<-TropicalPolytope.extrapolation.HAR_v4(B, x0, I = i,k=(d-1))
    proj<-project_pi(P,x)
    if(trop.dist(x,proj)<=1e-8){
      count<-count+1
      har_points1<-rbind(har_points1,x)
    }
    har_points[i,]<-x
    x0<-x
  }
  r<-count/sze
  VolB<-(d)*R^(d-1)
  VolP<-r*VolB
  return(list(r,VolB,VolP))
}

##########Drawing 3D hyperplanes for a given set of vertices ############
## Author :  David Barnhill                                            ##
## Date   :  November 21st 2022                                        ##
## Program:  This code produces the 3D tropical min-plus hyperplanes   ##
##           for a given set of vertices defining a tropical polytope. ##
## Input  :  A matrix, D, of vertices; maximum length, di; minimum     ##
##           value for each coordinate axis, mi; maximum value for each##
##           coordinate axis, ma; logical to determine whether to plot.##
## Output :  Plot of min-plus hyperplanes defining a tropical polytope.##
## Execute:  type in R as                                              ##  
#########################################################################
hypers3d_min<-function(D,di,mi,ma,plt=FALSE){
  cl<-rainbow(nrow(D))
  d<-dim(D)
  if((d[2]-1)==3){
    if(plt==TRUE){
      plot3d(D[,2],D[,3],D[,4],type='n',xlim=c(mi,ma),ylim=c(mi,ma),zlim=c(mi,ma),xlab='x2',ylab = 'x3',zlab='x4')
    }
    for (i in 1:(nrow(D))){
      x<-D[i,]
      xs1<-c(x[2],x[2]-di,x[2]-di,x[2])
      ys1<-c(x[3],x[3]-di,x[3]-di,x[3])
      zs1<-c(x[4],x[4]-di,x[4]+di,x[4]+di)
      
      xs2<-c(x[2],x[2]-di,x[2]-di,x[2])
      ys2<-c(x[3],x[3]-di,x[3]+di,x[3]+di)
      zs2<-c(x[4],x[4]-di,x[4]-di,x[4])
      
      xs3<-c(x[2],x[2]-di,x[2]+di,x[2]+di)
      ys3<-c(x[3],x[3]-di,x[3]-di,x[3])
      zs3<-c(x[4],x[4]-di,x[4]-di,x[4])
      
      xs4<-c(x[2],x[2]+di,x[2]+di,x[2])
      ys4<-c(x[3],x[3],x[3]+di,x[3]+di)
      zs4<-c(x[4],x[4],x[4],x[4])
      
      xs5<-c(x[2],x[2]+di,x[2]+di,x[2])
      ys5<-c(x[3],x[3],x[3],x[3])
      zs5<-c(x[4],x[4],x[4]+di,x[4]+di)
      
      xs6<-c(x[2],x[2],x[2],x[2])
      ys6<-c(x[3],x[3],x[3]+di,x[3]+di)
      zs6<-c(x[4],x[4]+di,x[4]+di,x[4])
      
      polygon3d(xs1,ys1,zs1,col=cl[i],coords = c(1,3))
      polygon3d(xs2,ys2,zs2,col=cl[i],coords=c(1,2))
      polygon3d(xs3,ys3,zs3,col=cl[i],coords=c(1,2))
      polygon3d(xs4,ys4,zs4,col=cl[i],coords = c(1,2),plot = TRUE)
      polygon3d(xs5,ys5,zs5,col=cl[i],coords = c(1,3),plot = TRUE)
      polygon3d(xs6,ys6,zs6,col=cl[i],coords = c(2,3),plot=TRUE)
    }
  }
  else{
    if(plt==TRUE){
      plot(D[,2],D[,3],type='n',xlim=c(mi,ma),ylim=c(mi,ma),xlab='x2',ylab = 'x3')
    }
    for (i in 1:(nrow(D))){
      x<-D[i,]
      xs1<-c(x[2],x[2]-di)
      ys1<-c(x[3],x[3]-di)
      
      xs2<-c(x[2],x[2]+di)
      ys2<-c(x[3],x[3])
      
      xs3<-c(x[2],x[2])
      ys3<-c(x[3],x[3]+di)
      
      lines(xs1,ys1,col=cl[i])
      lines(xs2,ys2,col=cl[i])
      lines(xs3,ys3,col=cl[i])
    }
  }
}

############Tropical Determinant of a matrix#############################
## Author :  David Barnhill                                            ##
## Date   :  January 8th 2023                                          ##
## Program:  This code calculates the tropical determinant of a square ##
##           matrix.  It assumes that the matrix is non-singular.      ##
## Input  :  A square matrix, P, of tropical coordinates with the      ##
##           the matrix arranged with the coordinates as row vectors.  ##
## Output(list) :  1) value of the tropical determinant; 2) matrix of  ##
##                 original coordinates reordered by value of          ##
##                 contribution to tropical determinant (small to      ##
##                 large).                                             ##
## Execute:  type in R as                                              ##  
#########################################################################
tdet<-function(P){
  dd<-dim(P)
  if(dd[[1]]!=dd[[2]]){
    print("Not a Square Matrix!")
  }
  else{
    PP<-matrix(0,nrow(P),0,TRUE)
    K<-t(P)
    i<-ncol(K)
    tdet<-0
    while (i>1){
      ind<-which.max(K[i,])
      tdet<-tdet+K[i,ind]
      PP<-cbind(K[,ind],PP)
      K<-K[,-ind]
      if(i>1){
        i=i-1
      }
    }
  
  tdet<-tdet+K[i]
  PP<-cbind(K,PP)%*%diag(1,nrow(PP))
  return(list(tdet,PP))}
}

#### Second Version ####
tdets<-function(P){
  dd<-dim(P)
  B<-P
  if(dd[[1]]!=dd[[2]]){
    return(print("Not a Square Matrix!"))
  }
  else{
    tds<-c(0,0)
    perms<-permn(dd[[1]])
    tdet<-0
    i=1
    max_ind<-0
    while(i<=length(perms)){
      k<-perms[[i]]
      t<-0
      for(j in 1:length(k)){
        t<-t+P[j,k[j]]
      }
      if(t>tdet){
        tdet<-t
        tds<-c(tdet,0)
        max_ind<-perms[[i]]
      }
      else if (t==tdet){
        tds<-c(tdet,t)
      }
      i=i+1
    }
    if(tds[1]==tds[2]){
      return(print('Singular'))
    }
    else{
      for (i in 1:length(max_ind)){
        P[max_ind[i],]<-B[i,]
        }
      return(list(tdet,P))
    }
  }
}


############Rounding Algorithm for a tropical polytope###################
## Author :  David Barnhill                                            ##
## Date   :  January 8th 2023                                          ##
## Program:  This code enumerates the pseudo-vertices of the boundary  ##
##           of a full-dimensional sector of a tropical polytope.      ##
## Input  :  A matrix, V, of vertices defining the tropical convex hull##
##           of a polytope P.                                          ##
## Output :  matrix of pseudo-vertices for the full dimensional part   ##
##           of a tropical polytope.                                   ##
## Execute:  type in R as                                              ##  
#########################################################################

rounding<-function(P){
  PP_star<-tdet(P)
  PP<-PP_star[[2]]
  for (i in 1:ncol(PP)){
    PP[,i]<-PP[,i]-PP[i,i]
  }
  b1<-rep(0,(nrow(PP)*(nrow(PP)-1)))
  k<-1
  for (i in 1:nrow(PP)){
    for(j in 1:ncol(PP)){
      if (i!=j){
        b1[k]<-PP[i,j]
        k=k+1}
    }
  }
  
  B<-matrix(NA,0,ncol(P))
  for (i in 1:nrow(PP)){
    p<-rep(0,ncol(PP))
    for (j in 1:ncol(PP)){
      pt<-p
      if(i==j) next
      if (j!=i){}
      pt[i]=1
      pt[j]=-1
      B<-rbind(B,pt)
    }
  }
  H<-makeH(-B[,-1],-b1)
  V<-scdd(H)
  V1<-cbind(rep(0,nrow(V$output[,c(3:ncol(V$output))])),V$output[,c(3:ncol(V$output))])
  return(V1)
}

#########Minimum Encompassing Ball for a tropical polytope###############
## Author :  David Barnhill                                            ##
## Date   :  January 9th 2023                                          ##
## Program:  This code produces a center point and radius of the       ##
##           minimum encompassing ball for a tropical polytope.        ##
## Input  :  A matrix, A, of vertices defining the tropical convex hull##
##           of a polytope P. Points are rows of matrix A.             ##
## Output :  Center point for a minimum encompassing ball for a        ##
##           tropical polytope, P, and associated radius.              ##
## Execute:  type in R as                                              ##  
#########################################################################

min_enc_ball<-function(A){
  P<-permutations(ncol(A),2)
  V<-matrix(0,0,ncol(A))
  bb<-c()
  for (j in 1:nrow(P)){
    for (i in 1:nrow(A)){
      k<-A[i,P[j,1]]-A[i,P[j,2]]
      bb<-append(bb,k)
      a<-rep(0,ncol(A))
      a[P[j,1]]<-1
      a[P[j,2]]<--1
      V<-rbind(V,a)
    }
  }
  r<-rep(-1,nrow(V))
  V<-cbind(V,r)
  f.obj<-c(rep(0,ncol(A)),1)
  f.con<-V
  f.dir <- rep("<=",nrow(V))
  f.rhs<-bb
  
  res<-lp ("min", f.obj, f.con, f.dir, f.rhs)
  sol<-res$solution
  cent<-sol[1:ncol(A)]
  rad<-sol[length(sol)]
  return(list(cent,rad))
}

#########Maximum Inscribed Ball for a tropical polytope##################
## Author :  David Barnhill                                            ##
## Date   :  January 10th 2023                                         ##
## Program:  This code produces a center point and radius of the       ##
##           maximum inscribed ball for a tropical polytope.           ##
## Input  :  A matrix, A, of vertices defining the tropical convex hull##
##           of a polytope P. Points are rows of matrix A.             ##
## Output :  Center point and radius for a minimum encompassing ball   ##
##           for a tropical polytope, P, and associated radius.        ##
## Execute:  type in R as                                              ##  
#########################################################################

max_ins_ball<-function(A){
  n<-rep(0,(ncol(A)))
  A<-tdets(A)[[2]]
  for (i in 2:ncol(A)){
    if(min(A[,i])<0){
      n[i]<-min(A[,i])
      A[,i]<-A[,i]-min(A[,i])
    }
  }
  W<-matrix(0,0,ncol(A))
  P<-permutations((ncol(A)),2)
  bb<-c()
  for(i in 1:nrow(A)){
    A[i,]<-A[i,]-A[i,i]
  }
  A<-t(A)
  for (i in 1:nrow(P)){
    if(P[i,1]==1){
      a<-rep(0,ncol(A))
      a[P[i,2]-1]=1
      a[length(a)]=2
      W<-rbind(W,a)
      bb<-append(bb,A[P[i,1],P[i,2]])
    }
    else if(P[i,2]==1){
      a<-rep(0,ncol(A))
      a[P[i,1]-1]<--1
      W<-rbind(W,a)
      bb<-append(bb,A[P[i,1],P[i,2]])
    }
    else{
      a<-rep(0,ncol(A))
      a[P[i,1]-1]<--1
      a[P[i,2]-1]<-1
      a[length(a)]<-1
      W<-rbind(W,a)
      bb<-append(bb,A[P[i,1],P[i,2]])
    }
  }
  f.con<-W
  f.dir<-c(rep("<=",(ncol(A)^2-ncol(A))))
  f.rhs<--bb
  f.obj<-c(rep(0,(ncol(A)-1)),1)
  sol<-lp ("max", f.obj, f.con, f.dir, f.rhs)
  solu<-sol$solution
  rad<-solu[length(solu)]
  cent<-solu[1:(length(solu)-1)]+n[2:length(n)]+rad
  return(list(rad,c(0,cent)))
}

