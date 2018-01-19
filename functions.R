# Number of neighbours (4) of pixel (a,b) with value x[i,j]=col
# 4 neighbours relation: (a-1,b),(a+1,b), (a,b-1) & (a,b+1)
xneig4=function(x,a,b,col){
  n=dim(x)[1];m=dim(x)[2]
  nei=c(x[a-1,b]==col,x[a,b-1]==col)
  if (a!=n)
    nei=c(nei,x[a+1,b]==col)
  if (b!=m) 
    nei=c(nei,x[a,b+1]==col)
  sum(nei) # True = 1, False=0
}

sumising=function(G=2,niter=10^3,n,m=n,beta){ 
  # Metropolis-Hastings Sampler
  # simulations by Gibbs with warmin stage
  S=0
  # Initialization:
  x=matrix(sample(c(1:G),n*m,replace=TRUE),n,m)
  # iteration:
  for (i in 1:niter){
    s=0
    # a random ordering of the elements of the lattice I
    sampl=sample(1:(n*m))
    for (k in 1:(n*m)){ # n*m = card(I)
      #print(c("iter:",i," px:",k))
      # current pixel: sampl[k] = px + (py-1)*n
        a = (sampl[k]-1)%%n+1 # (modulo) sampl[k] mod n = x-coord
        b = (sampl[k]-1)%/%n+1 # (integer division) sampl[k] \ n +1 = y-coord
        xcur=x[a,b] # idem x[sample[k]]
        # uniform proposal on the G-1 other possible values
        xtilde=sample((1:G)[-xcur],1) # uniform proposal 
        # log of the acceptance probability 
        lacpt = beta*(xneig4(x,a,b,xtilde) - xneig4(x,a,b,xcur))
        if (log(runif(1))< lacpt){
          x[a,b]=xtilde 
         }
        s=s+xneig4(x,a,b,x[a,b]) 
      }
    if (2*i>niter) #burning i > niter/2
      S=S+s
  }
  return(2*S/niter) # E(S(x))
}

# TRUNCATED NORMAL VARIATE, BASED ON INVERSE CDF:
truncnorm=function(n,mu,tau2,a,b)
{
  qnorm(pnorm(b,mu,sqrt(tau2))-runif(n)*(pnorm(b,mu,sqrt(tau2))-pnorm(a,mu,sqrt(tau2))),mu,sqrt(tau2))
}


# This function adresses the reconstruction of an image distributed from
# a Potts model based on a noisy version of this image. 
# The purpose of image segmentation (Chapter 8, Bayesian Essentials with R)
# is to cluster pixels into homogeneous classes without supervision or
# preliminary definition of those classes, based only on the spatial
# coherence of the structure. 
# The underlying algorithm is an hybrid Gibbs sampler.

reconstruct=function(G=6,niter,y) #number of Gibbs iterations; blurred image defined as a matrix y
{
  # y = log(ha/FUV) --> instead of grey levels [0..255]; [-1,3] for example
  numb=dim(y)[1]
  x=0*y
  mu=matrix(0,niter,G)
  sigma2=rep(0,niter)
  # init values, resorting to the background pixels (with values -10)
  mu[1,]= c(-10,quantile(y[which(y!=-10)],seq(0.01,0.99,length.out = G-1))) # c(0.05,0.25,0.50,0.75,0.95) 
  # In M77 there are 'outliers'
  for(k in 1:numb){
    for (l in 1:numb){
      if(y[k,l] > quantile(y[which(y!=-10)],0.99)) y[k,l] = quantile(y[which(y!=-10)],0.99)
    }
  }
  sigma2[1]= sd(y[which(y!=-10)])^2 #(sd(y)^2/G) #*10 #100
  beta=rep(1,niter)
  
  xcum=matrix(0,numb^2,G) # for each px of the image, account freq. of each color
  n=rep(0,G) # number of neighbors of each color 1<= g <= G(6)
  
  #Z=seq(0,2,length=21)
  #for (i in 1:21) Z[i]=sumising(niter=100,numb,beta=Z[i])
  thefunc=approxfun(seq(0,2,length=21),dali) # ,Z)
  # the vector dali comes from a crude if time
  # consumming approximation of the normalizing constant
  # based on 21 points
  i = 2
  #for (i in 2:niter)
  while(i<=niter & beta[i] < 2 & beta[i]>0)
  {
    
    lvr=0 # S(x)
    
    for (k in 1:numb)
    {
      for (l in 1:numb)
      {
        for(g in 1:G){
          n[g]=xneig4(x,k,l,g) #neighbors of px(k,l) with color g 
        }
        # full conditional distribution of x: 
        # p(x[k,l]=g |y,beta,sigma,mu) propto exp(beta*n - (1/2*sigma2)*(y[k,l]-mu[g])^2) 
        x[k,l]=sample(1:G,1,prob=exp(beta[i-1]*n)*dnorm(y[k,l],mu[i-1,],sqrt(sigma2[i-1])))
        #cat("probs =",exp(beta[i-1]*n)*dnorm(y[k,l],mu[i-1,],sqrt(sigma2[i-1])),"\n","k:",k," l:",l,"  sample = ",x[k,l],"\n")
        # account the no. of g coincidences for px(k,l) over the niters. (post distr of x)
        xcum[(k-1)*numb+l,x[k,l]]=xcum[(k-1)*numb+l,x[k,l]]+1
        # S(x) = # sum of neighbors with the same color, for whole lattice
        lvr=lvr+n[x[k,l]] 
      }
    }
    # full conditional distribution of mu: truncated normal distribution on [μg-1, μg+1] (setting μ0 = 0 and μG+1 = 255) 
    # with mean = sg/ng ; sg = sum of the observations y, allocated to category g; ng= number of observations of g
    # and variance = sigma^2/ng
    for(g in 1:G){
      if(g==1) col1 = min(y) else col1 = mu[i,g-1]
      if(g==G) col2 = max(y) else col2 = mu[i-1,g+1]
      if(empty(as.data.frame(y[x==g]))) mu[i,g]=runif(1,col1,col2)
      else  mu[i,g]=truncnorm(1,mean(y[x==g]),sqrt(sigma2[i-1]/sum(x==g)),col1,col2)
    }
    # full conditional distribution of sigma2: inv-gamma
    sese=sum((y-mu[i,1])^2*(x==1))
    for(g in 2:G)         
      sese = sese + sum((y-mu[i,g])^2*(x==g))
    sigma2[i]=1/rgamma(1,numb^2/2,sese/2)
    # MCMC sampler targeting the posterior distribution of beta:
    # proposal: uniform move with range 2h (h=0.05)
    repeat{
      betatilde=beta[i-1]+runif(1,-0.05,0.05)
      if(betatilde >0 & betatilde <2) break
    }
    print(c(i,betatilde))
    # acceptance ratio (log) associated with the pair (betatilde,beta)
    laccept=lvr*(betatilde-beta[i-1])+integrate(thefunc,betatilde,beta[i-1])$value
    if (runif(1)<=exp(laccept)) beta[i]=betatilde else beta[i]=beta[i-1]
    
    # print(round(c(i,beta[i],mu[i,],sigma2[i]),2))
    i = i+1
  }
  
  list(beta=beta,mu=mu,sigma2=sigma2,xcum=xcum)
}
