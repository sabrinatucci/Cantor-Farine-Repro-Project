# Figure S2: exclusivity surface plots for different initial connectivity and resource share type -------

## Load data on disproportionate share type: this scenario calculates individual payoffs in which the group-level resource share is determined by the size of the groups, considering that larger groups outcompete smaller ones. There are simulations for three initial network connectivity values: Tprob=0.2, Tprob=0.5, Tprob=0.8

load(paste(getwd(), "/data/1_exclusivity_surface.RData", sep=""))
source("setup.R")

# Preparing the data for surface plot
output <- output[which(!is.na(output[,4])),]
output <- as.data.frame(output)
output$ID <- 1:nrow(output)

# Bins for the plot
r.bins <- seq(min(output$r),max(output$r),1)
N.bins <- seq(min(output$N),max(output$N),1)

# figure labels
figlabel <- c("(A)", "(B)", "(C)")

par(mfcol=c(3,3), mar=c(1,1,1,1))

for(s in 1:3){
  
  # Data to plot: only cases when resource.size<= half of population
  input <- output[which(output$tprob==tprob[s] & output$r<=(output$N/2)),]
  
  # Fit surface
  surf <- locfit(Exclusivity~lp(N,r,nn=0.05,scale=F, h=0.1,deg=1), data=input)
  
  # Z-axis
  plotcol="black"
  zmax <- 1
  zmin <- 0
  z <- matrix(predict(surf,newdata=expand.grid(r=r.bins,N=N.bins),type="response"),nrow=length(r.bins),ncol=length(N.bins),byrow=FALSE)
  N.mat <- matrix(rep(N.bins,each=length(r.bins)),ncol=length(N.bins),nrow=length(r.bins))
  r.mat <- matrix(rep(r.bins,each=length(N.bins)),ncol=length(N.bins),nrow=length(r.bins),byrow=TRUE)
  z[which(r.mat > (N.mat/2))] <- NA
  z[which(z > 1)] <- 1
  z[which(z < 0)] <- 0
  minz <- min(z,na.rm=T)
  nrz <- nrow(z)
  ncz <- ncol(z)
  
  # Create colors
  nbcol <- 100
  jet.colors <- blue2green2red(nbcol)
  jet.colors2 <- add.alpha(jet.colors,alpha=0.6)
  # Compute the z-value at the facet centres
  zfacet <- (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
  # Recode facet z-values into color indices
  facetcol <- cut(zfacet,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
  zcol <- cut(z,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
  
  
  # Plot surface with transparent backpoints
  res <- persp(r.bins, N.bins, matrix(NA,nrow=length(r.bins),ncol=length(N.bins)), theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors[facetcol], ylab="",xlab="",zlab="",zlim=c(zmin,zmax),xlim=range(r.bins),ylim=range(N.bins),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75, main="")
  
  input2 <- input
  input2$r <- input2$r + rnorm(nrow(input),mean=0,sd=0.2)
  input2$N <- input2$N + rnorm(nrow(input),mean=0,sd=0.2)
  input2$pred <- predict(surf,newdata=input[,c(2,1)], type="response")
  
  input3 <- input2[which((input2$Exclusivity - input2$pred) < 0),]
  for (i in 1:nrow(input3)) {
    lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="darkgrey")
  }
  points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)
  par(new=TRUE)
  
  # add transparent surface
  res <- persp(r.bins, N.bins, z, theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors2[facetcol], ylab="",xlab="",zlab="Exclusivity",zlim=c(zmin,zmax),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75)
  input3 <- input2[which((input2$Exclusivity - input2$pred) >= 0),]
  for (i in 1:nrow(input3)) {
    lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="grey")
  }
  points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)
  #text(trans3d(45,5,1.2,res), paste(figlabel[s]," Tprob =", tprob[s]))
  text(trans3d(45,5,1.2,res), paste(figlabel[s]))

}







### Load data on equal share type: this scenario considers no disproportionate individual payoffs, i.e. there is equal allocation of resource to all individuals irrespective of their group size. There are simulations for three initial network connectivity values: Tprob=0.2, Tprob=0.5, Tprob=0.8

rm(list=ls())
load(paste(getwd(), "/data/1_exclusivity_surface_equalresourceshare.RData", sep=""))
source("setup.R")


# Preparing the data for surface plot
output <- output[which(!is.na(output[,4])),]
output <- as.data.frame(output)
output$ID <- 1:nrow(output)

# Bins for the plot
r.bins <- seq(min(output$r),max(output$r),1)
N.bins <- seq(min(output$N),max(output$N),1)

# figure labels
figlabel <- c("(D)", "(E)", "(F)")

for(s in 1:3){
  
  # Data to plot: only cases when resource.size<= half of population
  input <- output[which(output$tprob==tprob[s] & output$r<=(output$N/2)),]
  
  # Fit surface
  surf <- locfit(Exclusivity~lp(N,r,nn=0.05,scale=F, h=0.1,deg=1), data=input)
  
  # Z-axis
  plotcol="black"
  zmax <- 1
  zmin <- 0
  z <- matrix(predict(surf,newdata=expand.grid(r=r.bins,N=N.bins),type="response"),nrow=length(r.bins),ncol=length(N.bins),byrow=FALSE)
  N.mat <- matrix(rep(N.bins,each=length(r.bins)),ncol=length(N.bins),nrow=length(r.bins))
  r.mat <- matrix(rep(r.bins,each=length(N.bins)),ncol=length(N.bins),nrow=length(r.bins),byrow=TRUE)
  z[which(r.mat > (N.mat/2))] <- NA
  z[which(z > 1)] <- 1
  z[which(z < 0)] <- 0
  minz <- min(z,na.rm=T)
  nrz <- nrow(z)
  ncz <- ncol(z)
  
  # Create colors
  nbcol <- 100
  jet.colors <- blue2green2red(nbcol)
  jet.colors2 <- add.alpha(jet.colors,alpha=0.6)
  # Compute the z-value at the facet centres
  zfacet <- (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
  # Recode facet z-values into color indices
  facetcol <- cut(zfacet,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
  zcol <- cut(z,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
  
  
  # Plot surface with transparent backpoints
  res <- persp(r.bins, N.bins, matrix(NA,nrow=length(r.bins),ncol=length(N.bins)), theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors[facetcol], ylab="",xlab="",zlab="",zlim=c(zmin,zmax),xlim=range(r.bins),ylim=range(N.bins),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75, main="")
  
  input2 <- input
  input2$r <- input2$r + rnorm(nrow(input),mean=0,sd=0.2)
  input2$N <- input2$N + rnorm(nrow(input),mean=0,sd=0.2)
  input2$pred <- predict(surf,newdata=input[,c(2,1)], type="response")
  
  input3 <- input2[which((input2$Exclusivity - input2$pred) < 0),]
  for (i in 1:nrow(input3)) {
    lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="darkgrey")
  }
  points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)
  par(new=TRUE)
  
  # add transparent surface
  res <- persp(r.bins, N.bins, z, theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors2[facetcol], ylab="",xlab="",zlab="",zlim=c(zmin,zmax),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75)
  input3 <- input2[which((input2$Exclusivity - input2$pred) >= 0),]
  for (i in 1:nrow(input3)) {
    lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="grey")
  }
  points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)
  #text(trans3d(45,5,1.2,res), paste(figlabel[s]," Tprob =", tprob[s]))
  text(trans3d(45,5,1.2,res), paste(figlabel[s]))
}







#####
### When running the model for longer (n.reps=500) under equal allocation of resources the same pattern for disproportionate allocation of resources emerge ###
## the commented code below creates a restricted parameter space, runs the simulation with the edited model1() function and the edited simulation() function

## Load packages and functions
# source("setup.R")
#
## simplified parameter space, with fewer replicates 
# N = c( 20,  30,  40,  50,  70, 100, 140, 170, 200)
# r = c(5,  9, 13, 17, 25, 31, 41 ,51)
# tprob = c(0.2, 0.5, 0.8)
# replic = 100
# payoff.type = "equal"
# 
## parameter space
# variables <- expand.grid(N,r,tprob,payoff.type)
# colnames(variables) <- c("N","r","tprob","type")
# variables <- variables[which(variables$r<=(variables$N/2)),]      # only when resource.size <= half of population
# variables <- variables[rep(seq_len(nrow(variables)), each=replic),]
# variables <- split(variables, seq(nrow(variables)))
# 
### Loading edited model1 with coop.total after 250 t step
# model1 <- function(N, resource.size, n.reps, tprob, type="size-based", output="model"){
#   
#   # Setting initial directed network: start with a random set of edges; no loops
#   coop <- rgraph(N,tprob=tprob); diag(coop) <- NA
#   # Network with average number of foraging association across runs
#   coop.total <- matrix(0,nrow=N,ncol=N)  
#   # Mean individual payoffs
#   mean.payoff <- rep(NA,n.reps)
#   # All formed groups
#   groups.total <- list()
#   # Exclusivity
#   group.exclusivity <- rep(NA,n.reps)
#   
#   # Model run
#   for (zz in 1:n.reps) {
#     coop.NOW <- matrix(0,nrow=N,ncol=N)
#     
#     # extract pairs of cooperators (reciprocal links) out of the network
#     ABind <- cooperators(coop,N)
#     # form groups based on chain rule (A->B, B->C, then A->C)
#     groups <- identify.groups(ABind)
#     # calculate per capita payoff based on given resource patch size
#     if (length(groups) >= 1) {
#       payoffs <- calculate.payoffs(groups,resource.size,type) # calculate the group-level resource share based on group size
#       # individual payoffs
#       inds.payoffs <- rep(0,N)
#       for (i in 1:length(groups)) {
#         inds.payoffs[groups[[i]]] <- payoffs[i]
#         a <- expand.grid(groups[[i]],groups[[i]])
#         a <- a[which(a[,1] != a[,2]),]
#         coop.NOW[as.matrix(a)] <- 1  
#         if(zz>=250){
#           coop.total[as.matrix(a)] <- coop.total[as.matrix(a)]+1
#         }else{
#           coop.total[as.matrix(a)]
#         }
#       }
#     }
#     
#     #update network
#     coop <- update.coop.network(coop,inds.payoffs)
#     
#     if(output != "exclusivity"){
#       # Mean nonzero individual payoff
#       mean.payoff[zz] <- mean(inds.payoffs[which(inds.payoffs!=0)])
#       # Number and size of groups for each run
#       groups.total[[zz]] <- groups
#       # Exclusivity: proportion of times individuals were part of the foraging group for each run
#       group.exclusivity[zz] <- exclusivity(initial.size=N, runs=zz, patch=resource.size, all.ties=coop.total)
#     }
#   }
#   
#   # OUTPUTS
#   if(output != "exclusivity"){
#     group.number <- groups.number(all.groups=groups.total, runs=n.reps)
#     group.size <- average.group.size(groups.total)
#     if(output == "model"){
#       result.model <- list(coop.total, mean.payoff, group.number, group.size, group.exclusivity, coop.NOW, groups)
#       names(result.model) <- c("coop.total", "mean.payoffs", "n.groups", "mean.group.size", "exclusivity", "coop.NOW", "groups")
#       return(result.model)
#     } else { # output == 'sensitivity'
#       result.model <- cbind(1:n.reps, rep(N,n.reps), rep(resource.size,n.reps), rep(tprob,n.reps), group.exclusivity, mean.payoff, group.number, group.size)
#       colnames(result.model) <- c("time.step", "N","r","tprob","Exclusivity", "Mean.payoffs","n.groups", "Mean.group.size")
#       return(result.model)
#     }
#   } else { # output == "exclusivity"
#     group.exclusivity <- exclusivity(initial.size=N, runs=(n.reps-N), patch=resource.size, all.ties=coop.total)
#     result.model <- c(N, resource.size, tprob, group.exclusivity)
#     return(result.model)
#   }
# }
# 
### Loading simulation with n.reps=500
# simulation <- function(inputs, model, output){
#   cat("running parameters:", as.numeric(inputs),"\n")
# 
#   # assign input parameters
#   N <- as.numeric(inputs[1])
#   resource.size <- as.numeric(inputs[2])
#   tprob <- as.numeric(inputs[3])
#   type <- as.character(inputs[,4])
#   n.reps <- 500
# 
#   #  # when resource patch size is <= population size
#   #  if (resource.size <= N) { 
#   if(model=="model1") { 
#     result <- model1(N, resource.size, n.reps, tprob, type, output)
#   } else {
#     n.reps <- 1000
#     result <- model2(N, resource.size, n.reps, tprob, type, output)}
#   #}
#   return(result)
# }
#
#
### Running the simulation
# ptm <- proc.time()
# output <- do.call('rbind',lapply(variables, simulation, model="model1", output="exclusivity"))
# colnames(output) <- c("N","r","tprob","Exclusivity")
# cat(paste("simulation time:", round(((proc.time() - ptm)[3])/60, digits=2), "min"))
#
#####


# Instead: just load data on model1 exclusivity ran for longer (t=500), Tprob=0.2
rm(list=ls())
load(paste(getwd(), "/data/1_exclusivity_t250-tprob02-equal.RData", sep=""))
source("setup.R")

# figure labels
figlabel <- c("(G)", "(H)", "(I)")

# Plot surface
output <- output[which(!is.na(output[,4])),]
output <- as.data.frame(output)
output$ID <- 1:nrow(output)
r.bins <- seq(min(output$r),max(output$r),1)
N.bins <- seq(min(output$N),max(output$N),1)
s=1
{
    # Data to plot: only cases when resource.size<= half of population
    input <- output[which(output$tprob==tprob[s] & output$r<=(output$N/2)),]
    # Fit surface
    surf <- locfit(Exclusivity~lp(N,r,nn=0.09,scale=F, h=0.1,deg=1), data=input)
    # Z-axis
    plotcol="black"
    zmax <- 1
    zmin <- 0
    z <- matrix(predict(surf,newdata=expand.grid(r=r.bins,N=N.bins),type="response"),nrow=length(r.bins),ncol=length(N.bins),byrow=FALSE)
    N.mat <- matrix(rep(N.bins,each=length(r.bins)),ncol=length(N.bins),nrow=length(r.bins))
    r.mat <- matrix(rep(r.bins,each=length(N.bins)),ncol=length(N.bins),nrow=length(r.bins),byrow=TRUE)
    z[which(r.mat > (N.mat/2))] <- NA
    z[which(z > 1)] <- 1
    z[which(z < 0)] <- 0
    minz <- min(z,na.rm=T)
    nrz <- nrow(z)
    ncz <- ncol(z)
    # Create colors
    nbcol <- 100
    jet.colors <- blue2green2red(nbcol)
    jet.colors2 <- add.alpha(jet.colors,alpha=0.6)
    # Compute the z-value at the facet centres
    zfacet <- (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
    # Recode facet z-values into color indices
    facetcol <- cut(zfacet,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
    zcol <- cut(z,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
    res <- persp(r.bins, N.bins, matrix(NA,nrow=length(r.bins),ncol=length(N.bins)), theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors[facetcol], ylab="",xlab="",zlab="",zlim=c(zmin,zmax),xlim=range(r.bins),ylim=range(N.bins),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75, main="")
    input2 <- input
    input2$r <- input2$r + rnorm(nrow(input),mean=0,sd=0.2)
    input2$N <- input2$N + rnorm(nrow(input),mean=0,sd=0.2)
    input2$pred <- predict(surf,newdata=input[,c(2,1)], type="response")
    input3 <- input2[which((input2$Exclusivity - input2$pred) < 0),]
    for (i in 1:nrow(input3)) {
      lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="darkgrey")
    }
    points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)
    par(new=TRUE)
    res <- persp(r.bins, N.bins, z, theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors2[facetcol], ylab="",xlab="",zlab="",zlim=c(zmin,zmax),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75)
    input3 <- input2[which((input2$Exclusivity - input2$pred) >= 0),]
    for (i in 1:nrow(input3)) {
      lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="grey")
    }
    points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)
      text(trans3d(45,5,1.2,res), paste(figlabel[s]))

}




# Load data on model1 exclusivity ran for longer (t=500), Tprob=0.5
rm(list=ls())
load(paste(getwd(), "/data/1_exclusivity_t250-tprob05-equal.RData", sep=""))
source("setup.R")
figlabel <- c("(G)", "(H)", "(I)")

# Plot surface
output <- output[which(!is.na(output[,4])),]
output <- as.data.frame(output)
output$ID <- 1:nrow(output)
r.bins <- seq(min(output$r),max(output$r),1)
N.bins <- seq(min(output$N),max(output$N),1)
s=1
{
    # Data to plot: only cases when resource.size<= half of population
    input <- output[which(output$tprob==tprob[s] & output$r<=(output$N/2)),]
    # Fit surface
    surf <- locfit(Exclusivity~lp(N,r,nn=0.09,scale=F, h=0.1,deg=1), data=input)
    # Z-axis
    plotcol="black"
    zmax <- 1
    zmin <- 0
    z <- matrix(predict(surf,newdata=expand.grid(r=r.bins,N=N.bins),type="response"),nrow=length(r.bins),ncol=length(N.bins),byrow=FALSE)
    N.mat <- matrix(rep(N.bins,each=length(r.bins)),ncol=length(N.bins),nrow=length(r.bins))
    r.mat <- matrix(rep(r.bins,each=length(N.bins)),ncol=length(N.bins),nrow=length(r.bins),byrow=TRUE)
    z[which(r.mat > (N.mat/2))] <- NA
    z[which(z > 1)] <- 1
    z[which(z < 0)] <- 0
    minz <- min(z,na.rm=T)
    nrz <- nrow(z)
    ncz <- ncol(z)
    # Create colors
    nbcol <- 100
    jet.colors <- blue2green2red(nbcol)
    jet.colors2 <- add.alpha(jet.colors,alpha=0.6)
    # Compute the z-value at the facet centres
    zfacet <- (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
    # Recode facet z-values into color indices
    facetcol <- cut(zfacet,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
    zcol <- cut(z,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
    res <- persp(r.bins, N.bins, matrix(NA,nrow=length(r.bins),ncol=length(N.bins)), theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors[facetcol], ylab="",xlab="",zlab="",zlim=c(zmin,zmax),xlim=range(r.bins),ylim=range(N.bins),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75, main="")
    input2 <- input
    input2$r <- input2$r + rnorm(nrow(input),mean=0,sd=0.2)
    input2$N <- input2$N + rnorm(nrow(input),mean=0,sd=0.2)
    input2$pred <- predict(surf,newdata=input[,c(2,1)], type="response")
    input3 <- input2[which((input2$Exclusivity - input2$pred) < 0),]
    for (i in 1:nrow(input3)) {
      lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="darkgrey")
    }
    points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)
    par(new=TRUE)
    res <- persp(r.bins, N.bins, z, theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors2[facetcol], ylab="",xlab="",zlab="",zlim=c(zmin,zmax),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75)
    input3 <- input2[which((input2$Exclusivity - input2$pred) >= 0),]
    for (i in 1:nrow(input3)) {
      lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="grey")
    }
    points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)
      text(trans3d(45,5,1.2,res), paste(figlabel[2]))

}




# Load data on model1 exclusivity ran for longer (t=500), Tprob=0.8
rm(list=ls())
load(paste(getwd(), "/data/1_exclusivity_t250-tprob08-equal.RData", sep=""))
source("setup.R")
figlabel <- c("(G)", "(H)", "(I)")

# Plot surface
output <- output[which(!is.na(output[,4])),]
output <- as.data.frame(output)
output$ID <- 1:nrow(output)
r.bins <- seq(min(output$r),max(output$r),1)
N.bins <- seq(min(output$N),max(output$N),1)
s=1
{
    # Data to plot: only cases when resource.size<= half of population
    input <- output[which(output$tprob==tprob[s] & output$r<=(output$N/2)),]
    # Fit surface
    surf <- locfit(Exclusivity~lp(N,r,nn=0.09,scale=F, h=0.1,deg=1), data=input)
    # Z-axis
    plotcol="black"
    zmax <- 1
    zmin <- 0
    z <- matrix(predict(surf,newdata=expand.grid(r=r.bins,N=N.bins),type="response"),nrow=length(r.bins),ncol=length(N.bins),byrow=FALSE)
    N.mat <- matrix(rep(N.bins,each=length(r.bins)),ncol=length(N.bins),nrow=length(r.bins))
    r.mat <- matrix(rep(r.bins,each=length(N.bins)),ncol=length(N.bins),nrow=length(r.bins),byrow=TRUE)
    z[which(r.mat > (N.mat/2))] <- NA
    z[which(z > 1)] <- 1
    z[which(z < 0)] <- 0
    minz <- min(z,na.rm=T)
    nrz <- nrow(z)
    ncz <- ncol(z)
    # Create colors
    nbcol <- 100
    jet.colors <- blue2green2red(nbcol)
    jet.colors2 <- add.alpha(jet.colors,alpha=0.6)
    # Compute the z-value at the facet centres
    zfacet <- (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
    # Recode facet z-values into color indices
    facetcol <- cut(zfacet,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
    zcol <- cut(z,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
    res <- persp(r.bins, N.bins, matrix(NA,nrow=length(r.bins),ncol=length(N.bins)), theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors[facetcol], ylab="",xlab="",zlab="",zlim=c(zmin,zmax),xlim=range(r.bins),ylim=range(N.bins),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75, main="")
    input2 <- input
    input2$r <- input2$r + rnorm(nrow(input),mean=0,sd=0.2)
    input2$N <- input2$N + rnorm(nrow(input),mean=0,sd=0.2)
    input2$pred <- predict(surf,newdata=input[,c(2,1)], type="response")
    input3 <- input2[which((input2$Exclusivity - input2$pred) < 0),]
    for (i in 1:nrow(input3)) {
      lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="darkgrey")
    }
    points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)
    par(new=TRUE)
    res <- persp(r.bins, N.bins, z, theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors2[facetcol], ylab="Population size",xlab="Resource patch size",zlab="",zlim=c(zmin,zmax),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75)
    input3 <- input2[which((input2$Exclusivity - input2$pred) >= 0),]
    for (i in 1:nrow(input3)) {
      lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="grey")
    }
    points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)
      text(trans3d(45,5,1.2,res), paste(figlabel[3]))
}
