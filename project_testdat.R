# Test Data
dat <- read.csv("ClusterData/testdat.csv")

tProcess = function(){
  data = dat[1:15,1:2]
  # General PCA : to see the direction prescribed by general PCA
  M = cov(data)
  e = eigen(M)
  O = apply(e$vectors,2,function(x) x/sqrt(sum(x*x)))
  
  # Within Cluster Covariance Matrix
  W1 = cov(data[1:5,])  #cov of the first 5 points (cluster 1)
  W2 = cov(data[6:10,])
  W3 = cov(data[11:15,])
  W = (5-1)*(W1+W2+W3)/(15-1) #combined (the denominator is n-1)
  # In Between Cluster Covariance Matrix
  m1 = apply(data[1:5,],2,mean)    #Mean of the first cluster
  m2 = apply(data[6:10,],2,mean)
  m3 = apply(data[11:15,],2,mean)
  B = 5*cov(rbind(m1,m2,m3))   #Covariance matrix when all points in a
                               #cluster equals the cluster mean
  A = solve(W)%*%B # We needed to maximize X'BX/X'WX, which is the largest eig value of W^-1B i.e. A
  eig <- eigen(A)
  Q = apply(eig$vectors,2,function(x) x/sqrt(sum(x*x))) # Normalizing the eigen vectors
  scores = as.matrix(data) %*% as.matrix(Q) # Scores, i.e. data in the eigen basis
  scores = scores[,-2] # Last part of the data can be removed so as to support that the data is affine space of 1 dim
  return(list(centre=data,onb=Q,scores=scores,values=eig$values,non_clus_onb=O))
}

p = tProcess()
rangex = range(dat[1]) # Range of X values
rangey = range(dat[2]) # Range of y values
mean_data = apply(dat,2,mean) # Mean Data
plot(dat[1:5,],main="Comparing Regular PCA with PCA for Clusters",col="red",xlab="",ylab="",xlim = rangex,ylim = rangey) # First Cluster
par(new=TRUE) # So as to continue to plot on same window
plot(dat[6:10,],main="",col="blue",xlab="",ylab="",xlim = rangex,ylim = rangey) # Second Cluster
par(new=TRUE)
plot(dat[11:15,],main="",col="green",xlab="",ylab="",xlim = rangex,ylim = rangey) # Third Cluster
segments(x0=mean_data[1],y0=mean_data[2],x1=(mean_data[1]+p$non_clus_onb[1,1]),y1=(mean_data[2]+p$non_clus_onb[2,1]),col="blue",lwd=3) # direction of spread of data from Regular PCA
segments(x0=mean_data[1],y0=mean_data[2],x1=(mean_data[1]+p$onb[1,1]),y1=(mean_data[2]+p$onb[2,1]),col="red",lwd=3) # Modified direction of spread of data considering clusters
