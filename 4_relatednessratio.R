# Adding birth and death processes and keeping track of relatedness among agents

# Summary -----------------------------------------------------------------
#
# 1. Initial setup
# 2. Model run
# 3. Relatedness Ratio Surface







# 1. Initial setup --------------------------------------------------------

# Load packages and functions
source("setup.R")



# Define parameter space for surface plot
N <- c(20,30,40,50,60,70,80,90,100,120,140,160,180,200)  # Population/network size
r <- c(5,7,9,11,13,15,17,21,25,29,31,37,41,45,51)        # Resource patch size
tprob <- c(0.5)                                          # Initial network connectance
replic <- 100                                            # Number of model replicates
n.reps <- 1000                                           # Time steps
payoff.type <- 'size-based'                              # Type of group resource share

params.all <- expand.grid(N,r,tprob,payoff.type)
colnames(params.all) <- c("N","r","tprob","type")
params.all <- params.all[which(params.all$r<=(params.all$N/2)),]      # only when resource.size <= half of population
params.replic <- params.all[rep(seq_len(nrow(params.all)), each=replic),] #params.replic <- do.call(rbind, replicate(replic, params.all, simplify=FALSE))
variables <- split(params.replic, seq(nrow(params.replic)))



# Alternatively: load data from long simulations
#load(paste(getwd(), "/data/3_relatedness_surface.RData", sep=""))








# 2. Model run ------------------------------------------------------------


# Run model 2 across parameter space to calculate the relatedness ratio among members of foraging groups and non-members individuals at the end of the simulation. It  outputs a matrix with 5 columns ("N","r","tprob","Relatedness.ratio", logRelatedness.ratio) for the surface plot
ptm <- proc.time()
#output.surf <- do.call('rbind',lapply(variables, simulation, model="model2", output="sensitivity"))
output.surf <- matrix(NA, nrow(params.replic), 5)
for(i in 1:length(variables)){
  cat("parameters:", as.numeric(variables[[i]]), "; run",i,"out of ",length(variables),"\n")
  output.surf[i,] <- model2(N=as.numeric(variables[[i]][1]), resource.size=as.numeric(variables[[i]][2]), n.reps=n.reps, tprob=as.numeric(variables[[i]][3]), type=as.character(variables[[i]][,4]), output="sensitivity")
}
colnames(output.surf) <- c("N","r","tprob","Relatedness.ratio", "logRelatedness.ratio")
cat(paste("simulation time:", round(((proc.time() - ptm)[3])/60, digits=2), "min"))

# Saving data from simulation and all variables
save.image(file=paste(getwd(), "/data/3_relatedness_surface.RData", sep=""))







# 3. Relatedness ratio surface plots ---------------------------------------------



# 3.1. Log Relatedness ratio, Initial population size, initial connectance=0.5

# Data preparation
output.lrip <- output.surf
output.lrip[,1] <- params.replic[,1]
output.lrip <- output.lrip[which(!is.na(output.lrip[,5])),]
output.lrip <- output.lrip[which(is.finite(output.lrip[,5])),]
output.lrip <- as.data.frame(output.lrip)
output.lrip$ID <- 1:nrow(output.lrip)

input.lrip <- output.lrip 
#input,lrip <- output.lrip[which(output.lrip$tprob==tprob[s] & output.lrip$r<=(output.lrip$N/2)),]
#input.lrip <- output.lrip[which(output.lrip$tprob==tprob[s]),]

# surface
surf.lrip <- locfit(logRelatedness.ratio~lp(N,r,nn=0.05,scale=F, h=0.1,deg=1), data=input.lrip) 

# plot
plotsurf(input=input.lrip, surf=surf.lrip, yaxis="Initial population size",zaxis="Log Relatedness Ratio",xaxis="Resource patch size")




# 3.2. Relatedness ratio, Final population size, initial connectance=0.5

# Data preparation
output.rfp <- output.surf
output.rfp <- output.rfp[which(!is.na(output.rfp[,4])),]
output.rfp <- output.rfp[which(is.finite(output.rfp[,4])),]
output.rfp <- as.data.frame(output.rfp)
output.rfp$ID <- 1:nrow(output.rfp)

input.rfp <- output.rfp 
#input.rfp <- output.rfp[which(output.rfp$r<=(output.rfp$N/2)),]
#input.rfp <- output.rfp[which(output.rfp$tprob==tprob[s]),]

# surface
surf.rfp <- locfit(Relatedness.ratio~lp(N,r,nn=0.05,scale=F, h=0.1,deg=1), data=input.rfp) 

# plot
plotsurf(input=input.rfp, surf=surf.rfp, yaxis="Final population size",zaxis="Relatedness Ratio",xaxis="Resource patch size")