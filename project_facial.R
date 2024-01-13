# We start with a bunch of facial photographs, 5 photos of 10 persons each. 
# All the images are 200$\times$180 in size. So they are 36000 dimensional 
# vectors. Our aim is to reduce the dimension as much as possible without 
# losing our ability to recognise the faces.

# We start with loading the images.

library(jpeg) # Required library

loadImages = function() { # The function loads all the images
  name = c("male.pacole","male.pspliu","male.sjbeck","male.skumar","male.rsanti",
           "female.anpage","female.asamma","female.klclar","female.ekavaz","female.drbost")
    
  x = matrix(0,nrow=10*5,ncol=200*180) # X is the data matrix
  k = 0 # K is the counter 
  
  for(i in 1:10) { # i represents which person's photograph is considered 
    for(j in 1:5) { # j represents what numbered photograph is considered for the fixed person (by i)
      k = k +  # hence k is the overall counter
      #The first argument below should be the location of the face folder.
      filename = paste("E:/ISI/Sem 2/VM Project/grayfaces/",name[i],".",j,".jpg",sep="")
      x[k,] = as.vector(readJPEG(filename)) # so the photo is converted into a vector of dim 36,000
    }
  }
  return(x)
}

# Next, we perform PCA. The following function does precisely this, but without using R's built-in tools like
# prcomp or princomp, we perform the operations explicitly so that you can appreciate what goes on
# "behind the scene".

process = function(x) {
  meanx = apply(x,2,mean) # We are getting the mean picture
  y = scale(x,scale=F) # Rows of y are the cases, origin shifted to mean
  # Now as origin is shifted to mean, so second order raw moment is variance, hence the matrix Y'Y is Covariance matrix of x
  # But next we need eigen values of the same, and calculating eig values of 36000x36000 matrix is impractical.
  # Hence we use the trick that AB and BA have the nonzero eig values in common, and YY' is only 50x50 instead of 36000x36000,
  # so we compute YY' and take the eig values of the same.
  A = y %*% t(y) 
  eig = eigen(A) # Eig vector associated with largest eig value maximize x'Ax, and thats what we use in PCA
  
  P = t(y) %*% eig$vec[,-50] # Now, columns of P are eig vectors of Y'Y (Y'YY'x = lambda*Y'x), last vector is not that important as the total is just an affine space in 49 dimension 
  Q = apply(P,2,function(x) x/sqrt(sum(x*x))) #Columns of Q form onb for rowspace of Y, just scaling the orthogonal basis to onb
  
  scores = y %*% Q # scores[i,] is the tuple representation of i-th image in the eigen basis
  # score_cluster = apply(scores,1,function(x) Re(sum(x))) # A suitable norm to represent the picture
  # plot(score_cluster) # Plotting the norm
  # We see that the clusters are not clear, as we have not spent any affort to minimise the distance between clusters,
  # So next we are going to do that.
  
  # Now comes the part of the in-between and within covariance matrices 
  # To start, we only consider the first 10 eigen-vectors : claiming that
  # 10 dimension is a good estimate for the dimension spanned by the facial data
  
  # Hence we have the truncated data:
  trunc_onb = Q[,1:10] # Only first 10 eig vectors (ordered in decreasing value of corresponding eig value)
  trunc_data = y %*% trunc_onb
  
  # Within Covariance:
  W = (5-1)*(cov(trunc_data[1:5,])+cov(trunc_data[6:10,])+cov(trunc_data[11:15,])+cov(trunc_data[16:20,])+cov(trunc_data[21:25,])+cov(trunc_data[26:30,])+cov(trunc_data[31:35,])+cov(trunc_data[36:40,])+cov(trunc_data[41:45,])+cov(trunc_data[46:50,]))/(50-1)
  # In between Covariance:
  B = 5*cov(rbind(apply(trunc_data[1:5,],2,mean),apply(trunc_data[6:10,],2,mean),apply(trunc_data[11:15,],2,mean),apply(trunc_data[16:20,],2,mean),apply(trunc_data[21:25,],2,mean),apply(trunc_data[26:30,],2,mean),apply(trunc_data[31:35,],2,mean),apply(trunc_data[36:40,],2,mean),apply(trunc_data[41:45,],2,mean),apply(trunc_data[46:50,],2,mean)))
  
  # We want to maximize x'Bx/x'Wx, which is max eig value of W'B, hence :
  M = solve(W)%*%B # Solve(W) is W^-1
  new_eig <- eigen(M) # Eigen value and vectors of M
  new_coord = apply(new_eig$vectors,2,function(x) x/sqrt(sum(x*x))) # Again normalizing to get the onb 
  new_scores = trunc_data %*% new_coord # New scores are new tuple representations of the pictures in the new onb
  # print(new_scores)
  new_score_clusters = apply(new_scores,1,function(x) Re(sum(x))) # A suitable norm to identify the cluster of vectors (10 dim)
  plot(new_score_clusters) # Plotting the norm 
  
  
  return(list(centre = meanx, onb = Q, scores = scores, values = eig$values, new_onb = trunc_onb, new_scores = new_scores, new_norm = new_score_clusters, new_values = new_eig$values))
}

# Each principal component is again a 200*180 dimensional vector, and hence may
# be considered as an image. It might be instructive to take a look at these.
# The following function helps you to just that.

showFace = function(newCoord, i) {
  plot(1:2,ty='n',main="0")
  y = abs(newCoord$onb[,i])
  extreme = range(y)
  y = (y-extreme[1])/(extreme[2]-extreme[1])
  dim(y) = c(200,180)
  rasterImage(as.raster(y),1,1,2,2)
}

# The next function reconstucts the faces from the principal components.
# It starts with the "mean face" and then gradually adds the details.
# You'll need to hit "enter" to step through the process.

showSteps = function(newCoord,i) {
  meanx = newCoord$centre
  Q = newCoord$onb
  scores = newCoord$scores
  values = newCoord$values
  expl = 100*cumsum(newCoord$values)/sum(newCoord$values)
  
  coeff = as.vector(scores[i,])
  
  plot(1:2,ty='n',main="0")
  y = meanx
  dim(y) = c(200,180)
  rasterImage(as.raster(y),1,1,2,2)
  readline()
  for(k in 1:49) {
    if(k==1)
      temp = Q[,1]*coeff[1]
    else
      temp=Q[,1:k] %*% as.vector(coeff[1:k])
    recons = meanx + temp
    recons[recons<0]=0
    recons[recons>1]=1
    dim(recons) = c(200,180)
    plot(1:2,ty='n',main=paste(k,":",values[k],", ",expl[k]))
    rasterImage(as.raster(recons),1,1,2,2)
    readline()
  }
}


x = loadImages()
newcoord = process(x)
showSteps(newcoord,3)

