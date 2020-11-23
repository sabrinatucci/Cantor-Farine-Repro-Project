# Figure 1: surface plot for emerget properties of model 1 (number of groups, mean group size, individual payoff) and sample slice histograms networks across the representative parameter space

# Figure 1 -------------------------------------------------------------

# Load data and functions
load(paste(getwd(), "/data/2_sensitivity_simulation.RData", sep=""))
source("setup.R")

# Plotting Fig1 for # N=40, r=5, tprob=0.2


# generating the probability data and plotting surfaces
z40.5n = z.plot(input.data=output.sensi, variable="n.groups", plot.parameters=c(40,5,0.2), x.limit=40, y.limits=c(1,100))
z40.5g = z.plot(input.data=output.sensi, variable="group.size", plot.parameters=c(40,5,0.2), x.limit=40, y.limits=c(1,100))
z40.5p = z.plot(input.data=output.sensi, variable="payoffs", plot.parameters=c(40,5,0.2), x.limit=NA, y.limits=c(1,100))

# Plotting as 3D histograms (for groups up to 20 members, since a group is defined by >1 individual and N=40)
par(mfrow=c(3,1), mar=c(1,1,1,1))

res1 <- hist3D(x=1:20, y=1:100, z=z40.5n[1:20,], theta=120, phi=30, shade=0.2,zlab="Probability",ylab='Time',xlab='Number of groups', scale=T, add=F, colkey=F,ticktype = "detailed", las=1, alpha=0.3,nticks=5)
text(trans3d(20,0,1,res1), "(a)")

res2 <- hist3D(x=1:20, y=1:100,z=z40.5g[1:20,], theta=120, phi=30, shade=0.2,zlab="Probability",ylab='Time',xlab='Average group size', scale=T, add=F, colkey=F,ticktype = "detailed", las=1, alpha=0.3,nticks=5)
text(trans3d(20,0,1,res2), "(b)")

res3 <- hist3D(x=seq(from=0, to=5, by=0.6),y=1:100, z=z40.5p[1:9,], theta=120, phi=30, shade=0.2,zlab="Probability",ylab='Time',xlab='Average individual payoff', scale=T, add=F, colkey=F,ticktype = "detailed", las=1, alpha=0.3,nticks=5)
text(trans3d(1,0,1.1,res3), "(c)")




# Choosing time steps:
time=c(5,20,100) 

# Number of groups
par(mfrow=c(length(time)*3,1), mar=c(4,4,1,1))
for(i in time){ 
  plot(0:40, z40.5n[,i], type='l',xlab="", ylab="", las=1, ylim=c(0,1),xaxt='n', yaxt='n')
  axis(side = 1, at = seq(0, 40, by=10), labels = FALSE)
  axis(side = 2, at = seq(0.2, 1, by=0.4), labels = TRUE, las=1)
  text(37,0.9, paste('T =', i))
  if(i == time[2]){ title(ylab="Probability")}
  if(i == time[3]){ title(xlab="Number of groups")
                    axis(side = 1, at = seq(0, 40, by=10), labels = TRUE)}    
}

for(i in time){ 
  plot(0:40, z40.5g[,i], type='l',xlab="", ylab="", las=1, ylim=c(0,1),xaxt='n', yaxt='n')
  axis(side = 1, at = seq(0, 40, by=10), labels = FALSE)
  axis(side = 2, at = seq(0.2, 1, by=0.4), labels = TRUE, las=1)
  text(37,0.9, paste('T =', i))
  if(i == time[2]){ title(ylab="Probability")}
  if(i == time[3]){ title(xlab="Average group size")
                    axis(side = 1, at = seq(0, 40, by=10), labels = TRUE)}    
}

for(i in time){ 
  plot(seq(0,(dim(z40.5p)[1]-1)/2,0.5), z40.5p[,i], type='l',xlab="", ylab="", las=1, ylim=c(0,1),xaxt='n', yaxt='n')
  axis(side = 1, at = seq(0,(dim(z40.5p)[1]-1)/2,1), labels = FALSE)
  axis(side = 2, at = seq(0.2, 1, by=0.4), labels = TRUE, las=1)
  text(3.5,0.9, paste('T =', i))
  #vertical line to show the centralisation on payoffs=1.
  lines(c(1,1),c(0,1),lty=2,lwd=0.8)
  if(i == time[2]){ title(ylab="Probability")}
  if(i == time[3]){ title(xlab="Average individual payoff")
                    axis(side = 1, at = seq(0,(dim(z40.5p)[1]-1)/2,1), labels = TRUE)}    
}
