# Figure S4: Relatedness ratio surface plots variations---------------------------------------------

# Load data and functions
load(paste(getwd(), "/data/3_relatedness_surface.RData", sep=""))
source("setup.R")

par(mfrow=c(1,2), mar=c(1,2,1,1))

s=1

# Fig S4 (a) Relatedness ratio, Final population size, initial connectance=0.5

# Data preparation
output.rfp <- output.surf
#output.rfp[,1] <- params.replic[,1] #initial population size
output.rfp <- output.rfp[which(!is.na(output.rfp[,4])),]
output.rfp <- output.rfp[which(is.finite(output.rfp[,4])),]
output.rfp <- as.data.frame(output.rfp)
output.rfp$ID <- 1:nrow(output.rfp)

input.rfp <- output.rfp 
#input.rfp <- output.rfp[which(output.rfp$r<=(output.rfp$N/2)),]
#input.rfp <- output.rfp[which(output.rfp$tprob==tprob[s]),]

# surface
surf.rfp <- locfit(Relatedness.ratio~lp(N,r,nn=0.2,scale=F, h=0.1,deg=1), data=input.rfp) 

# plot
res1 <- plotsurf(input=input.rfp, surf=surf.rfp, yaxis="Final population size",zaxis="Relatedness Ratio",xaxis="Resource patch size")
text(trans3d(48,5,0.45,res1), "(a)")



# FigS4 (b) Log Relatedness ratio, Final population size, initial connectance=0.5

# Data preparation
output.lrip <- output.surf
#output.lrip[,1] <- params.replic[,1] # initial population size
output.lrip <- output.lrip[which(!is.na(output.lrip[,5])),]
output.lrip <- output.lrip[which(is.finite(output.lrip[,5])),]
output.lrip <- as.data.frame(output.lrip)
output.lrip$ID <- 1:nrow(output.lrip)

input.lrip <- output.lrip 
#input.lrip <- output.lrip[which(output.lrip$tprob==tprob[s] & output.lrip$r<=(output.lrip$N/2)),]
#input.lrip <- output.lrip[which(output.lrip$tprob==tprob[s]),]

# surface
surf.lrip <- locfit(logRelatedness.ratio~lp(N,r,nn=0.20,scale=F, h=0.1,deg=1), data=input.lrip) 

# plot
res2 <- plotsurf(input=input.lrip, surf=surf.lrip, yaxis="Final population size",zaxis="Log Relatedness Ratio",xaxis="Resource patch size")
text(trans3d(48,5,3,res2), "(b)")
