library(ggplot2)
library(plyr)
source('funtions.R')

Z=seq(0,2,by=0.1) # discretize beta range [0,2]
# Constant of Potts model, image of size 300, G=6 homog. colors/clusters
for (i in 1:length(Z)) Z[i]=sumising(G=6,niter=100,n=300,beta=Z[i])
par(mfrow=c(1,1))
plot(seq(0,2,by=0.1),Z)
dali = Z

# PREPARE DATA as INPUT:
# Reading data
# 3 cols (px,py,flux)
galha <- as.data.frame(read.table('data/m83_ha_extR.list')) # Halpha obs. flux
galfuv <- as.data.frame(read.table('data/m83_fuv_extR.list')) # FUV obs. flux

#the size of the images is 370x370, backgrund = -10
gal = data.frame(V1=rep(1:370,times=370),V2=rep(1:370,each=370),V3=-10)
ind = which(galha$V3> 0 & galfuv$V3> 0)
gal[ind,]$V3 = log10(galha[ind,]$V3) -log10(galfuv[ind,]$V3) # log(Halpha/FUV)

#plotting the flux ratio map
ggplot(data=gal[which(gal$V3>-10),])+geom_raster(aes(V1, V2,fill=V3))

# Take a sub-image of size 300x300, same as dali's size
ind = which(gal$V1 %in% (1:300) & gal$V2 %in% (51:350)) #M83
galcut = gal[ind,]
# ggplot(data=galcut[which(galcut$V3 > -10),])+geom_raster(aes(V1, V2,fill=V3))

# give data the corresponding format as the input image for reconstruction funtion
m83 = matrix(galcut$V3,300,300)

# IMAGE SEGMENTATION into G=6 (5 homogeneous regions + 1 the background )
titus=reconstruct(G=6,niter=20,m83)

#returns a permutation which rearranges its first argument into ascending order
# keep the largest: mode or max of the post. distr. of x
affect=function(u) { order(u)[6] }

aff=apply(titus$xcum,1,affect)
aff=t(matrix(aff,300,300))
#Plotting the result
image(1:300,1:300,aff,col=gray(5:1/5),xlab="",ylab="")


