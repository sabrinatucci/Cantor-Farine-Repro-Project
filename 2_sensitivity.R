# Sensitivity analysis for model 1 parameters: how model outputs (number of foraging groups, average group size, average individual payoff) vary with model inputs (network size, resource patch size, connectance)

# Summary -----------------------------------------------------------------
#
# 1. Initial setup
# 2. Model run
# 3. Results: 3D plots
# 4. Results: Slice histograms
# 5. Saving working space




# 1. Initial setup --------------------------------------------------------

# Load packages and functions
source("setup.R")

# Define parameter space
N <- c(40,120)                  # Population/network size
r <- c(5,15,25,35)              # Resource patch size
tprob <- c(0.2)                 # Initial network connectance
replic <- 500                   # Number of model replicates
payoff.type <-'size-based'      # Type of group resource share

variables <- expand.grid(N,r,tprob,payoff.type)
variables <- variables[rep(seq_len(nrow(variables)), each=replic),]
variables <- split(variables, seq(nrow(variables)))




# 2. Model run ------------------------------------------------------------

# Run model1 across the parameter space, and outputs matrix with 8 columns("time.step", "N","r","tprob","Exclusivity", "Mean.payoffs","n.groups", "Mean.group.size")
ptm <- proc.time()
output.sensi <- do.call('rbind',lapply(variables, simulation, model="model1",output="sensitivity"))
# Run time (min)
cat(paste("Total runs:", length(variables), "; Simulation time:", round(((proc.time() - ptm)[3])/60, digits=2), "min"))


# Alternatively: load simulation data
#load(paste(getwd(), "/data/2_sensitivity_simulation.RData", sep=""))




# 3. Results: 3D plots --------------------------------------

# Plotting probability of Number of groups, Average Group size, and Average individual payoff versus total simulation time, for different sets of parameters

par(mfrow=c(4,3), mar=c(1,1,1,1))

# N=40, r=5, tprob=0.2
z40.5n = z.plot(input.data=output.sensi, variable="n.groups", plot.parameters=c(40,5,0.2), x.limit=10, y.limits=c(1,100))
z40.5g = z.plot(input.data=output.sensi, variable="group.size", plot.parameters=c(40,5,0.2), x.limit=20, y.limits=c(1,100))
z40.5p = z.plot(input.data=output.sensi, variable="payoffs", plot.parameters=c(40,5,0.2), x.limit=NA, y.limits=c(1,100))

# N=40, r=15, tprob=0.2
z40.15n = z.plot(input.data=output.sensi, variable="n.groups", plot.parameters=c(40,15,0.2), x.limit=10, y.limits=c(1,100))
z40.15g = z.plot(input.data=output.sensi, variable="group.size", plot.parameters=c(40,15,0.2), x.limit=20, y.limits=c(1,100))
z40.15p = z.plot(input.data=output.sensi, variable="payoffs", plot.parameters=c(40,15,0.2), x.limit=10, y.limits=c(1,100))

# N=120, r=25, tprob=0.2
z120.25n = z.plot(input.data=output.sensi, variable="n.groups", plot.parameters=c(120,25,0.2), x.limit=10, y.limits=c(1,100))
z120.25g = z.plot(input.data=output.sensi, variable="group.size", plot.parameters=c(120,25,0.2), x.limit=20, y.limits=c(1,100))
z120.25p = z.plot(input.data=output.sensi, variable="payoffs", plot.parameters=c(120,25,0.2), x.limit=10, y.limits=c(1,100))

# N=120, r=35, tprob=0.2
z120.35n = z.plot(input.data=output.sensi, variable="n.groups", plot.parameters=c(120,35,0.2), x.limit=10, y.limits=c(1,100))
z120.35g = z.plot(input.data=output.sensi, variable="group.size", plot.parameters=c(120,35,0.2), x.limit=20, y.limits=c(1,100))
z120.35p = z.plot(input.data=output.sensi, variable="payoffs", plot.parameters=c(120,35,0.2), x.limit=10, y.limits=c(1,100))







# 4. Slice histograms -----------------------------------------------------


# Plotting propability profile of Number of groups, Average Group size, and Average individual payoff in a given time step

# Choosing time steps:5,10,20,50,100
time=c(5,10,20,50,100) 

# Number of groups
# N=40, r=5, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(0:40, z40.5n[,i], type='l',xlab="Number of groups", ylab="", main=paste('Time step =', i), las=1, ylim=c(0,1))
  if(i == time[1]){ title(ylab="Probability")}    
}

# N=40, r=15, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(0:40, z40.15n[,i], type='l',xlab="Number of groups", ylab="", main=paste('Time step=', i), las=1, ylim=c(0,1))
  if(i == time[1]){ title(ylab="Probability")}    
}

# N=120, r=25, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(0:120, z120.25n[,i], type='l',xlab="Number of groups", ylab="", main=paste('Time step=', i), las=1, ylim=c(0,1))
  if(i == time[1]){ title(ylab="Probability")}    
}

# N=120, r=35, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(0:120, z120.35n[,i], type='l',xlab="Number of groups", ylab="", main=paste('Time step=', i), las=1, ylim=c(0,1))
  if(i == time[1]){ title(ylab="Probability")}    
}



# Average group size
# N=40, r=5, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(0:40, z40.5g[,i], type='l',xlab="Average group size", ylab="", main=paste('Time step=', i), las=1, ylim=c(0,1))
  if(i == time[1]){ title(ylab="Probability")}    
}

# N=40, r=15, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(0:40, z40.15g[,i], type='l',xlab="Average group size", ylab="", main=paste('Time step=', i), las=1, ylim=c(0,1))
  if(i == time[1]){ title(ylab="Probability")}    
}

# N=120, r=25, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(0:120, z120.25g[,i], type='l',xlab="Average group size", ylab="", main=paste('Time step=', i), las=1, ylim=c(0,1))
  if(i == time[1]){ title(ylab="Probability")}    
}

# N=120, r=35, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(0:120, z120.35g[,i], type='l',xlab="Average group size", ylab="", main=paste('Time step=', i), las=1, ylim=c(0,1))
  if(i == time[1]){ title(ylab="Probability")}    
}


# Individual payoffs
# N=40, r=5, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(seq(0,(dim(z40.5p)[1]-1)/2,0.5), z40.5p[,i], type='l',xlab="Average individual payoff", main=paste('Time step=', i), ylab="", las=1, ylim=c(0,1))
  #vertical line to show the centralisation on payoffs=1.
  lines(c(1,1),c(0,1),lty=2,lwd=0.8)
  if(i == time[1]){ title(ylab="Probability")}    
}


# N=40, r=15, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(seq(0,(dim(z40.15p)[1]-1)/2,0.5), z40.15p[,i], type='l',xlab="Average individual payoff", main=paste('Time step=', i), ylab="", las=1, ylim=c(0,1))
  #vertical line to show the centralisation on payoffs=1.
  lines(c(1,1),c(0,1),lty=2,lwd=0.8)
  if(i == time[1]){ title(ylab="Probability")}    
}

# N=120, r=25, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(seq(0,(dim(z120.25p)[1]-1)/2,0.5), z120.25p[,i], type='l',xlab="Average individual payoff", main=paste('Time step=', i), ylab="", las=1, ylim=c(0,1))
  #vertical line to show the centralisation on payoffs=1.
  lines(c(1,1),c(0,1),lty=2,lwd=0.8)
  if(i == time[1]){ title(ylab="Probability")}    
}

# N=120, r=35, tprob=0.2
par(mfrow=c(1,length(time)), mar=c(4,4,1,1))
for(i in time){ 
  plot(seq(0,(dim(z120.35p)[1]-1)/2,0.5), z120.35p[,i], type='l',xlab="Average individual payoff", main=paste('Time step=', i), ylab="", las=1, ylim=c(0,1))
  #vertical line to show the centralisation on payoffs=1.
  lines(c(1,1),c(0,1),lty=2,lwd=0.8)
  if(i == time[1]){ title(ylab="Probability")}    
}




# 5. 3D slice histograms (example plots)

time2=c(5,15,25,35,45,55,65,75,85,95)

auxn=matrix(0,nrow(z40.5n), ncol(z40.5n)); auxn[,time2]=z40.5n[,time2]
auxg=matrix(0,nrow(z40.5g), ncol(z40.5g)); auxg[,time2]=z40.5g[,time2]
auxp=matrix(0,nrow(z40.5p), ncol(z40.5p)); auxp[,time2]=z40.5p[,time2]

par(mfrow=c(3,1), mar=c(1,1,1,1))

hist3D(z=auxn[1:20,], space = c(0,0.9), theta=120, phi=30, shade=0.2,zlab="probability",ylab='time',xlab='number of groups')
hist3D(z=auxg[1:20,], space = c(0,0.9), theta=120, phi=30, shade=0.2,zlab="probability",ylab='time',xlab='average group size')
hist3D(z=auxp, space = c(0,0.9), theta=120, phi=30, shade=0.2,zlab="probability",ylab='time',xlab='average individual payoff')




# 5. Saving working space  ---------------------------------------------------------------

# saving data from simulation and all variables
save.image(file=paste(getwd(), "/data/2_sensitivity_simulation.RData", sep=""))
