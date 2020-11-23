# Figure 4: Relatedness ratio surface plot and relatedness from different areas of the parameter space



# Figure 4 ---------------------------------------------

# Load data and functions
load(paste(getwd(), "/data/3_relatedness_surface.RData", sep=""))
source("setup.R")


# Fig4 (a): Log Relatedness ratio, Initial population size, initial connectance=0.5

par(mfrow=c(1,1), mar=c(1,2,1,1))

# Data preparation
output.lrip <- output.surf
output.lrip[,1] <- params.replic[,1]
output.lrip <- output.lrip[which(!is.na(output.lrip[,5])),]
output.lrip <- output.lrip[which(is.finite(output.lrip[,5])),]
output.lrip <- as.data.frame(output.lrip)
output.lrip$ID <- 1:nrow(output.lrip)

input.lrip <- output.lrip 

# surface
surf.lrip <- locfit(logRelatedness.ratio~lp(N,r,nn=0.2,scale=F, h=0.1,deg=1), data=input.lrip) 

# plot
plotsurf(input=input.lrip, surf=surf.lrip, yaxis="Initial population size",zaxis="Log Relatedness Ratio",xaxis="Resource patch size")





# Fig 4 (b-e): Plotting relatedness rates, from  100 model replicates

# load data and functions
load(paste(getwd(), "/data/3_relatedness_Fig4.RData", sep=""))
source("setup.R")

### Alternatively, run the simulation
## Define representative parameter space
#Nr <- c(60,160)    # Population/network size
#rr <- c(7,31)      # Resource patch size
#tprobr <- c(0.5)   # Initial network connectance
#n.reps <- 1000     # Time steps
#replicr <- 100     # Number of model replicates
#
#params <- expand.grid(Nr,rr,tprobr)
#params.r <- params[rep(seq_len(nrow(params)), each=replicr),]
#params.r <- split(params, seq(nrow(params)))
#
## running the model
#output.relat.r <- list()
#output.relat.r <- lapply(params.r, simulation, model="model2",output="relatedness")
#
## Organizing model replicates; every 100 replicates of a given set of params: output.relat.replic[[params]][[replic]][[time.step]][[output]]
#output.relat.replic.r <- list(output.relat.r[c(1:100)], output.relat.r[c(101:200)], output.relat.r[c(201:300)], output.relat.r[c(301:400)])
###



# Plot relatedness among Foraging Group Members, non-members and random across representative parameter space
par(mfrow=c(2,2), mar=c(4,4,2,0.5))
lim=c(0.2,0.2,0.1,0.1)
for(j in 1:nrow(params)){
  # relatedness y-axis limits
  lim=numeric()
  for(z in 1:length(output.relat.replic.r[[j]])){
    aux=output.relat.replic.r[[j]][[z]]
    lim[z]=max(aux[[5]]$all.relatedness[,1:3], na.rm=T)
  }
  # empty plot
  plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(0, round(max(lim),digits=2)), las=1, main=paste("N=", params[j,1], ", r=", params[j,2], ', c=', params[j,3], ', replic=100',sep=''))
#  plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(0, lim[j]), las=1, main=paste("N=", params[j,1], ", r=", params[j,2], ', c=', params[j,3], ', replic=100',sep=''))
  # overlay members, non-members, random relatedness
  for(i in 1:length(output.relat.replic.r[[j]])){
    # output.relat.replic.r[[params]][[replic]][[time.step]][[output]]
    points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.random.relatedness, pch=20, col=rgb(red=1,green=0,blue=0,alpha=0.01))
    points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.noncoop.relatedness, pch=20, col=rgb(red=0,green=0,blue=1,alpha=0.02))
    points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.members.relatedness, pch=20, col=rgb(red=0,green=0,blue=0,alpha=0.3))
  }
  if(j == 2) legend('topright', c("Random", "Non-members", "Members"), cex=0.9, col=c("red","blue", "black"), bty="n", lwd=c(3, 3)) 
  if(j %in% c(1,3)){ title(ylab="Mean Relatedness")}
  if(j %in% c(3,4)){ title(xlab="Time step")}   
}




# Plot Log relatedness ratio across representative parameter space

# First: Calculate mean and 95%CI of log relatedness ratio per time step, across model replicates
# pulling out relatedness ratio of each time step, from all 100 model replicats (rows=time steps; cols=ratio from model replicates) 
replicratio <- rep(list(matrix(NA, nrow=n.reps, ncol=replicr)), nrow(params))
for(j in 1:nrow(params)){
  for(i in 1:replicr){
    replicratio[[j]][,i] = output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$relat.ratio
  }
}
# calculate logratio mean and 95%CI (rows=time steps; cols=mean(logratio), 2.5% & 97.5%)
ratiostat <- rep(list(matrix(NA, nrow=n.reps, ncol=3)), nrow(params))
for(j in 1:nrow(params)){
  for(t in 1:n.reps){
    aux=log(replicratio[[j]][t,])
    aux=aux[is.finite(aux)]
    aux=aux[!is.na(aux)]
    ratiostat[[j]][t,1] <- mean(aux)
    ratiostat[[j]][t,2:3] <- quantile(aux, probs=c(0.025, 0.975), type=2)
  }
}

par(mfrow=c(2,2), mar=c(4,4,2,0.5))
for(j in 1:nrow(params)){
  # relatedness y-axis limits
  ulim <- llim <- numeric()
  for(z in 1:length(output.relat.replic.r[[j]])){
    aux=output.relat.replic.r[[j]][[z]]
    lim=log(aux[[5]]$all.relatedness$relat.ratio)
    lim[!is.finite(lim)] <- NA
    ulim[z]=max(lim, na.rm=T, is.finite=T)
    llim[z]=min(lim, na.rm=T, is.finite=T)     
  }
  # empty plot
  plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(round(min(llim,na.rm=T),digits=1), round(max(ulim,na.rm=T),digits=1)), las=1, main=paste("N=", params[j,1], ", r=", params[j,2], ', c=', params[j,3], ', replic=100',sep=''))
#  plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(-4,6), las=1, main=paste("N=", params[j,1], ", r=", params[j,2], ', c=', params[j,3], ', replic=100',sep=''))
  # overlay relatedness ratios
  for(i in 1:length(output.relat.replic.r[[j]])){
    points(1:n.reps, log(output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$relat.ratio), pch=20, col=rgb(red=0,green=0,blue=0,alpha=1))
  }
    x <- roll_mean(1:n.reps,n=20)
    y <- roll_mean(ratiostat[[j]][,1],n=20,na.rm=TRUE)
    y1 <- roll_mean(ratiostat[[j]][,2],n=20,na.rm=TRUE)
    y2 <- roll_mean(ratiostat[[j]][,3],n=20,na.rm=TRUE)
    a <- predict(loess(y ~ x, span=0.05),newdata=data.frame(x))
    a1 <- predict(loess(y1 ~ x, span=0.05),newdata=data.frame(x))
    a2 <- predict(loess(y2 ~ x, span=0.05),newdata=data.frame(x))
    # overlay 95%CI polygon and mean line
    polygon(c(x,rev(x)),
            c(a1,rev(a2)),
            col=rgb(red=0,green=1,blue=0.1,alpha=0.5),border=NA)
    lines(x,a, col=rgb(red=0.1,green=1,blue=0.1,alpha=1))
    abline(0,0, col=rgb(red=0,green=0,blue=0,alpha=1), lty=2)

  if(j %in% c(1,3)){ title(ylab="Log Mean Relatedness Members/Random")}
  if(j %in% c(3,4)){ title(xlab="Time step")}   
}






## Plot each part of the Figure 4 separately to be edited in Illustrator


#(b)
par(mfrow=c(1,2), mar=c(0.5,3,0.5,0.5), las=1)
lim=c(0.2,0.2,0.2,0.2)
j=1
params[j,]
#plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(0, lim[j]), las=1)
# remove all axis labels to be edited in Adobe
plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(0, lim[j]), las=1, xaxt='n', yaxt='n')
axis(side = 1, at = seq(0, n.reps, by=200), labels = FALSE)
axis(side = 2, at = seq(0, lim[j], by=0.05), labels = FALSE)
# overlay members, non-members, random relatedness
for(i in 1:length(output.relat.replic.r[[j]])){
  # output.relat.replic.r[[params]][[replic]][[time.step]][[output]]
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.random.relatedness, pch=20,col=rgb(red=1,green=0,blue=0,alpha=0.03))
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.noncoop.relatedness, pch=20,col=rgb(red=0,green=0,blue=1,alpha=0.03))
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.coop.relatedness, pch=20,col=rgb(red=0,green=0,blue=0,alpha=0.5))
}
# overlay relatedness ratios
#plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(-4,6), las=1)
plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(-4,6), las=1, xaxt='n', yaxt='n')
axis(side = 1, at = seq(0, n.reps, by=200), labels = FALSE)
axis(side = 2, at = seq(-4, 6, by=2), labels = FALSE)
for(i in 1:length(output.relat.replic.r[[j]])){
  points(1:n.reps, log(output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$relat.ratio), pch=20,col=rgb(red=0,green=0,blue=0,alpha=0.5))
}
x <- roll_mean(1:n.reps,n=20)
y <- roll_mean(ratiostat[[j]][,1],n=20,na.rm=TRUE)
y1 <- roll_mean(ratiostat[[j]][,2],n=20,na.rm=TRUE)
y2 <- roll_mean(ratiostat[[j]][,3],n=20,na.rm=TRUE)
a <- predict(loess(y ~ x, span=0.05),newdata=data.frame(x))
a1 <- predict(loess(y1 ~ x, span=0.05),newdata=data.frame(x))
a2 <- predict(loess(y2 ~ x, span=0.05),newdata=data.frame(x))
# overlay 95%CI polygon and mean line
polygon(c(x,rev(x)),
        c(a1,rev(a2)),
        col=rgb(red=0,green=1,blue=0.1,alpha=0.5),border=NA)
lines(x,a, col=rgb(red=0.1,green=1,blue=0.1,alpha=1))
abline(0,0, col=rgb(red=0,green=0,blue=0,alpha=1), lty=2)






#(c)
par(mfrow=c(2,1), mar=c(0.5,3,0.5,0.5), las=1)
j=2
params[j,]
#plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(0, lim[j]), las=1)
# remove all axis labels to be edited in Adobe
plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(0, lim[j]), las=1, xaxt='n', yaxt='n')
axis(side = 1, at = seq(0, n.reps, by=200), labels = FALSE)
axis(side = 2, at = seq(0, lim[j], by=0.05), labels = FALSE)
# overlay members, non-members, andom relatedness
for(i in 1:length(output.relat.replic.r[[j]])){
  # output.relat.replic.r[[params]][[replic]][[time.step]][[output]]
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.random.relatedness, pch=20,col=rgb(red=1,green=0,blue=0,alpha=0.03))
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.noncoop.relatedness, pch=20,col=rgb(red=0,green=0,blue=1,alpha=0.03))
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.coop.relatedness, pch=20,col=rgb(red=0,green=0,blue=0,alpha=0.5))
}
# overlay relatedness ratios
#plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(-4,6), las=1)
plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(-4,6), las=1, xaxt='n', yaxt='n')
axis(side = 1, at = seq(0, n.reps, by=200), labels = FALSE)
axis(side = 2, at = seq(-4, 6, by=2), labels = FALSE)
for(i in 1:length(output.relat.replic.r[[j]])){
  points(1:n.reps, log(output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$relat.ratio), pch=20,col=rgb(red=0,green=0,blue=0,alpha=0.5))
}
x <- roll_mean(1:n.reps,n=20)
y <- roll_mean(ratiostat[[j]][,1],n=20,na.rm=TRUE)
y1 <- roll_mean(ratiostat[[j]][,2],n=20,na.rm=TRUE)
y2 <- roll_mean(ratiostat[[j]][,3],n=20,na.rm=TRUE)
a <- predict(loess(y ~ x, span=0.05),newdata=data.frame(x))
a1 <- predict(loess(y1 ~ x, span=0.05),newdata=data.frame(x))
a2 <- predict(loess(y2 ~ x, span=0.05),newdata=data.frame(x))
# overlay 95%CI polygon and mean line
polygon(c(x,rev(x)),
        c(a1,rev(a2)),
        col=rgb(red=0,green=1,blue=0.1,alpha=0.5),border=NA)
lines(x,a, col=rgb(red=0.1,green=1,blue=0.1,alpha=1))
abline(0,0, col=rgb(red=0,green=0,blue=0,alpha=1), lty=2)






#(d)
par(mfrow=c(1,2), mar=c(0.5,3,0.5,0.5), las=1)
j=4
params[j,]
#plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(0, lim[j]), las=1)
# remove all axis labels to be edited in Adobe
plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(0, lim[j]), las=1, xaxt='n', yaxt='n')
axis(side = 1, at = seq(0, n.reps, by=200), labels = FALSE)
axis(side = 2, at = seq(0, lim[j], by=0.05), labels = FALSE)
# overlay members, non-members, random relatedness
for(i in 1:length(output.relat.replic.r[[j]])){
  # output.relat.replic.r[[params]][[replic]][[time.step]][[output]]
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.random.relatedness, pch=20,col=rgb(red=1,green=0,blue=0,alpha=0.03))
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.noncoop.relatedness, pch=20,col=rgb(red=0,green=0,blue=1,alpha=0.03))
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.coop.relatedness, pch=20,col=rgb(red=0,green=0,blue=0,alpha=0.5))
}
# overlay relatedness ratios
#plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(-4,6), las=1)
plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(-4,6), las=1, xaxt='n', yaxt='n')
axis(side = 1, at = seq(0, n.reps, by=200), labels = FALSE)
axis(side = 2, at = seq(-4, 6, by=2), labels = FALSE)
for(i in 1:length(output.relat.replic.r[[j]])){
  points(1:n.reps, log(output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$relat.ratio), pch=20,col=rgb(red=0,green=0,blue=0,alpha=0.5))
}
x <- roll_mean(1:n.reps,n=20)
y <- roll_mean(ratiostat[[j]][,1],n=20,na.rm=TRUE)
y1 <- roll_mean(ratiostat[[j]][,2],n=20,na.rm=TRUE)
y2 <- roll_mean(ratiostat[[j]][,3],n=20,na.rm=TRUE)
a <- predict(loess(y ~ x, span=0.05),newdata=data.frame(x))
a1 <- predict(loess(y1 ~ x, span=0.05),newdata=data.frame(x))
a2 <- predict(loess(y2 ~ x, span=0.05),newdata=data.frame(x))
# overlay 95%CI polygon and mean line
polygon(c(x,rev(x)),
        c(a1,rev(a2)),
        col=rgb(red=0,green=1,blue=0.1,alpha=0.5),border=NA)
lines(x,a, col=rgb(red=0.1,green=1,blue=0.1,alpha=1))
abline(0,0, col=rgb(red=0,green=0,blue=0,alpha=1), lty=2)





#(e)
par(mfrow=c(2,1), mar=c(0.5,3,0.5,0.5), las=1)
j=3
params[j,]
#plot(1, type="n", xlab="", ylab="Mean Relatedness", xlim=c(0, n.reps), ylim=c(0, lim[j]), las=1)
# remove all axis labels to be edited in Adobe
plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(0, lim[j]), las=1, xaxt='n', yaxt='n')
axis(side = 1, at = seq(0, n.reps, by=200), labels = FALSE)
axis(side = 2, at = seq(0, lim[j], by=0.05), labels = FALSE)
# overlay members, non-members, random relatedness
for(i in 1:length(output.relat.replic.r[[j]])){
  # output.relat.replic.r[[params]][[replic]][[time.step]][[output]]
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.random.relatedness, pch=20,col=rgb(red=1,green=0,blue=0,alpha=0.03))
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.noncoop.relatedness, pch=20,col=rgb(red=0,green=0,blue=1,alpha=0.03))
  points(1:n.reps, output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$mean.coop.relatedness, pch=20,col=rgb(red=0,green=0,blue=0,alpha=0.5))
}
# overlay relatedness ratios
#plot(1, type="n", xlab="Time step", ylab="Log Mean Relatedness Members/Random", xlim=c(0, n.reps), ylim=c(-4,6), las=1)
plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(-4,6), las=1, xaxt='n', yaxt='n')
axis(side = 1, at = seq(0, n.reps, by=200), labels = FALSE)
axis(side = 2, at = seq(-4, 6, by=2), labels = FALSE)
for(i in 1:length(output.relat.replic.r[[j]])){
  points(1:n.reps, log(output.relat.replic.r[[j]][[i]][[5]]$all.relatedness$relat.ratio), pch=20,col=rgb(red=0,green=0,blue=0,alpha=0.5))
}
x <- roll_mean(1:n.reps,n=20)
y <- roll_mean(ratiostat[[j]][,1],n=20,na.rm=TRUE)
y1 <- roll_mean(ratiostat[[j]][,2],n=20,na.rm=TRUE)
y2 <- roll_mean(ratiostat[[j]][,3],n=20,na.rm=TRUE)
a <- predict(loess(y ~ x, span=0.05),newdata=data.frame(x))
a1 <- predict(loess(y1 ~ x, span=0.05),newdata=data.frame(x))
a2 <- predict(loess(y2 ~ x, span=0.05),newdata=data.frame(x))
# overlay 95%CI polygon and mean line
polygon(c(x,rev(x)),
        c(a1,rev(a2)),
        col=rgb(red=0,green=1,blue=0.1,alpha=0.5),border=NA)
lines(x,a, col=rgb(red=0.1,green=1,blue=0.1,alpha=1))
abline(0,0, col=rgb(red=0,green=0,blue=0,alpha=1), lty=2)
