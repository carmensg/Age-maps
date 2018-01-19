# Libraries required
library(parallel)
require(monmlp)
require(plyr)

# Loading the perceptron ----------------------------------------------
load("perceptron.Rdata") # loading perceptron

# Functions needed to do posterior samples ------------------------------
prior <- function(){
  z <- runif(1,0.001,0.04)
  fqh <- runif(1,0.1,1) 
  age <- 19.9*rbeta(n = 1, shape1 = 1, shape2 = 1)+0.1 # escalada 0.1-20
  while(monmlp.predict(t(as.matrix(c(z,fqh,age))), ann)[3] >=1){
    #imf <- sample(1:3, size = 1)
    z <- runif(1,0.001,0.04)
    fqh <- runif(1,0.1,1) 
    age <- 19.9*rbeta(n = 1, shape1 = 1, shape2 = 1)+0.1 # escalada 0.1-20
  }
  theta <-c(z,fqh,age) #c(imf, z, age)
  return(theta)
}

loglikelihood <- function(theta, ratio.obs){#,sigma1,sigma2,rho){
  mu <- monmlp.predict(t(as.matrix(theta)), ann)
  mu1 <- 10^mu[1] #  Ha.true
  sigma1 <- sqrt(0.05^2+0.04^2)*mu1 # Assume a relative error of 5%+4%
  mu2 <- 10^mu[2] # Fuv.true: flux density!!
  sigma2 <- sqrt(0.25^2+0.04^2)*mu2 # Assume a relative error of 25%+4%
  rho <- mu[3]
  z <- ratio.obs
  
  a <- sqrt( z^2/sigma1^2-2*rho*z/(sigma1*sigma2) + 1/sigma2^2 )
  b <- mu1*z/sigma1^2 - rho*(mu1+mu2*z)/(sigma1*sigma2) + mu2/sigma2^2
  c <- mu1^2/sigma1^2 - 2*rho*mu1*mu2/(sigma1*sigma2) + mu2^2/sigma2^2
  d <- exp((b^2-c*a^2)/(2*(1-rho^2)*a^2))
  q <- b/(sqrt(1-rho^2)*a)
  
  lik <- b*d/(sqrt(2*pi)*sigma1*sigma2*a^3)*( 2*pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)-1)+
    sqrt(1-rho^2)/(pi*sigma1*sigma2*a^2)* exp(-c/(2*(1-rho^2)))
  return(log(lik))
}

# Functions needed to do Nested Sampling --------------------------------

NestedSampling <- function(ratio.obs, prior. = prior, loglikelihood. = loglikelihood,
                           dim.par = 3, M = 100, N = 100){
  # SET PRIOR OBJECTS
  S <- matrix(NA, N, dim.par) 
  for(j in 1:N){
    S[j,]  <- prior.()
  }
  loglik <- adply(S,.margins = 1, loglikelihood., ratio.obs)
  S <- as.data.frame(S) 
  # Añado una 4º col. "loglik" a S, que es la col. V1 de loglik: Li
  S$loglik <- loglik$V1
  # Ahora S contiene los pares (theta_i) y L_i 
  
  # SETTING UP THE logZ AND H, THEN nested.object :
  H <- 0 # Information, initially 0
  logZ <- -.Machine$double.xmax # ln(Evidence Z, initially 0)
  # nested.object es una lista, con 3 elementos: S, logZ y H  --> nested.object$H por ejem.
  nested.object <- list(S = S, logZ = logZ, H = H) # SURVIVORS
  inactive.points <- data.frame() #data frame with 0 columns and 0 rows : DISCARDED
  
  # SETTING UP LOOP TERMINATION
  i <- 1
  cond1 <- i < M
  f <- -2 # As Farhan Feroz in doi:10.1111/j.1365-2966.2007.12353.x (0.2 in log-evidence) ! LOG10 NO LN !!!!!
  max.loglik <- max(nested.object$S$loglik) # dentro del elem. S del nested.object, tomo max de su 4º col: L_i 
  logXj <- 0
  # If max(L1,...L_N)*Xj < f*Zj --> termination
  cond2 <- (max.loglik + logXj) > f + logZ  # "f" realmente es logf
  
  # NESTED SAMPLING LOOP: 
  while(cond1 & cond2) {
    # Evolve copied object within constraints
    evolved.nested.object <- NestedSamplingStep(nested.object, i, prior.,
                                                loglikelihood., ratio.obs)
    nested.object <- evolved.nested.object$no
    inactive.points <- rbind(inactive.points, evolved.nested.object$ip)   
    # setting up termination conditions
    max.loglik <- max(nested.object$S$loglik)
    logXj <- -i/N # standard crude assignment (shrink interval)
    logZ <- nested.object$logZ    
    #print(paste("Iteration = ", i, "worstloglik = ", evolved.nested.object$ip$loglik,"logZ=", logZ,sep = " "))
    cond2 <- (max.loglik + logXj) > f + logZ
    i <- i + 1
    cond1 <- i < M
  }
  # Ordering by likelohood the active points or nested objtects, S:
  nested.object$S <- nested.object$S[order(nested.object$S$loglik),]
  # Complete Z:
  # Simply:
  nested.object$logZ <- log.plus(nested.object$logZ,log(sum(exp(nested.object$S$loglik))/N*exp(logXj)))
  NS <- rbind(inactive.points,nested.object$S)
  NS$nest <- 1:dim(NS)[1]
  NS$post <- exp(NS$loglik-NS$nest/N-nested.object$logZ)
  return(list(NS=NS,logZ=nested.object$logZ,H=nested.object$H))
  #return(list(no = nested.object, ip=inactive.points))
}

# Worst object in collection, with Weight=width*Likelihood
NestedSamplingStep <- function(nested.object, iteration, prior.,
                               loglikelihood., ratio.obs) {
  
  S <- nested.object$S
  H <- nested.object$H
  logZ <- nested.object$logZ
  # Record the lowest of the current likelihood values:
  worstloglik <- min(S$loglik)
  worstpoint <- which.min(S$loglik)  
  
  inactive.point <- S[worstpoint,]
  
  # Updating nested object
  #print(paste("Iteration = ", iteration, "loglik = ", worstloglik,sep = " "))
  N <- dim(S)[1]
  logwidth <- log(0.5*(exp(-(iteration-1)/N)-exp(-(iteration+1)/N))) # trapezoidal rule: log(Wi) = log((X(i-1)-X(i+1))/2)
  logWt <- logwidth + worstloglik # Wt = Li*Wi
  # Update evidence and information:
  logZnew <- log.plus(logZ, logWt) #logZnew <- log(Z+Wt)
  if(exp(worstloglik)>0){ # For cases where loglik = -infty !!
    H <- exp(logWt - logZnew) * worstloglik + exp(logZ - logZnew) * (H + logZ) - logZnew
    #print(paste("Iteration = ", iteration, "loglik = ", worstloglik,"H = ",H,"logZ= ",logZ,sep = " "))
  }
  # Kill worst object in favour of copy of different survivor: 
  new.active.point <- SamplingNewCandidate(worstloglik,
    prior., loglikelihood., ratio.obs)
  S[worstpoint,] <- new.active.point
  nested.object$S <- S
  nested.object$logZ <- logZnew
  nested.object$H <- H
  
  #print(paste("Number Iterations = ", iteration,sep = " "))
  return(list(no = nested.object, ip = inactive.point))
}


SamplingNewCandidate <- function(worstloglik, prior., loglikelihood., ratio.obs){
    theta <- prior.()
    loglik <- loglikelihood.(theta, ratio.obs)
    
    if(loglik > worstloglik){
      return(c(theta, loglik))
    }else{
      cond <- TRUE
      while(cond){
        theta <- prior.()
        loglik <- loglikelihood.(theta, ratio.obs)
        cond <- loglik < worstloglik
      }
      return(c(theta, loglik))
    }
}

log.plus <- function(x,y) {
  if(x>y) x+log(1+exp(y-x))
  else    y+log(1+exp(x-y))
}

# Loading dataset -------------------------------------------------------

m83 <- as.matrix(read.table('m83.list')) # 370*370 = 136900 px
colnames(m83) <- c("px","py","ratio","Ha","Fuv","age1","age2") 
index <- which(m83[,4] >0 & m83[,5] >0)
ratio.obs <- as.numeric(10^(m83[index, "ratio"]))

# Doing the analysis ----------------------------------------------------
post.galaxy <- parLapply(cl, ratio.obs, NestedSampling)
save(post.galaxy, file = "m83post.Rdata")

# Doing the posterior analysis ----------------------------------------------------
m83$ModeAge <- rep(c(0),times=dim(m83)[1])
m83$MedianAge <- rep(c(0),times=dim(m83)[1])
m83$P25Age <- rep(c(0),times=dim(m83)[1])
m83$P75Age <- rep(c(0),times=dim(m83)[1])
m83$Zmode <- rep(c(0),times=dim(m83)[1])
m83$IMFmode <- rep(c(0),times=dim(m83)[1])

for(m in 1:length(index)){
  post.pixel = as.data.frame(post.galaxy[m])#post.galaxy[[indextoplay[n]]]
  colnames(post.pixel) <- c("Z","fqh","age","loglik","nest","post","logZ","H")
  sample <- post.pixel[unique(sample(nrow(post.pixel), size = 1000,replace = TRUE, prob = post.pixel$post)),]
  d <- density(sample$age)
  m83$ModeAge[total.pixels[m]] <- d$x[which.max(d$y)] # Mode
  m83$P25Age[total.pixels[m]] <- quantile(sample$age)[2] #P25
  m83$MedianAge[total.pixels[m]] <- quantile(sample$age)[3] #P50
  m83$P75Age[total.pixels[m]] <- quantile(sample$age)[4] #P75
  d <- density(sample$Z)
  m83$Zmode[total.pixels[m]] <- d$x[which.max(d$y)] # Mode Z
  d <- density(sample$fqh)
  m83$IMFmode[total.pixels[m]] <- d$x[which.max(d$y)] # Mode IMF
}

write.table(gal,"m83_postagemaps.txt")


