# Libraries required
require(monmlp)
require(plyr)

# Loading the perceptron ----------------------------------------------
load("perceptron.Rdata") # loading perceptron

# Loading Nested Sampling function ------------------------------
source("NestedSampling.R")

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


