# Adding birth and death processes and keeping track of relatedness among agents

# Summary -----------------------------------------------------------------
#
# 1. Initial setup
# 2. Model run
# 3. Relatedness Ratio
# 4. Relatedness Distribution
# 5. Network plots
# 6. Number of generations






# 1. Initial setup --------------------------------------------------------

# Load packages and functions
source("setup.R")


# Define representative parameter space
N <- c(40,120)                 # Population/network size
r <- c(5,15,25,35)             # Resource patch size
tprob <- c(0.5)                # Initial network connectance
replic <- 100                  # Number of model replicates
payoff.type <-'size-based'     # Type of group resource share



params <- expand.grid(N,r,tprob,payoff.type)
colnames(params) <- c("N","r","tprob","type")
#params <- params[which(params$r<=(params$N/2)),]      # only when resource.size <= half of population
params.replic.r <- params[rep(seq_len(nrow(params)), each=replic),]
variables.r <- split(params.replic.r, seq(nrow(params.replic.r)))

# Define representative parameter space
#Nr <- c(40,160)    # Population/network size
#rr <- c(7,31)      # Resource patch size
#tprobr <- c(0.5)  # Initial network connectance
#replicr <- 100           # Number of model replicates
#payoff.typer <-'size-based'     # Type of group resource share

#params <- expand.grid(Nr,rr,tprobr,payoff.typer)
#params.r <- params.r[rep(seq_len(nrow(params.r)), each=replicr),]
#params.r <- split(params.r, seq(nrow(params.r)))

# Representative areas of the parameter space for plotting
#params <- params.r[c(1,13,4,16,17,29,20,32,33,45,36,48),]



# Alternatively: load data from long simulations
#load(paste(getwd(), "/data/3_relatedness_samples.RData", sep=""))








# 2. Model run ------------------------------------------------------------


# 2.1 Typical model 2 runs

# Outputs final time step
network2f = model2(N=40,resource.size=5,n.reps=1000,tprob=0.2,type="size-based",output="model")
# Outputs for time steps 200, 400, 600, 800, 1000
network2s = model2(N=40,resource.size=5,n.reps=1000,tprob=0.2,type="size-based",output="samples")
# Outputs for time steps 200, 400, 600, 800, 1000, including all relatedness measures (slower)
network2 = model2(N=40,resource.size=5,n.reps=1000,tprob=0.2,type="size-based",output="relatedness")
# Outputs only relatedness ratio for the last time step for sensitivity analysis
network2r = model2(N=40,resource.size=5,n.reps=1000,tprob=0.2,type="size-based",output="sensitivity")




# 2.2 Run model2 with relatedness across the representative parameter space. For each replicate, it output the same list with 10 objects, for 5 steps along the simulaion time : coop.NOW, relatedness.link, relatedness.network, groups, relatedness.link.total, relatedness.network.total, ids, alive matrix, coop.total, all.relatedness
output.relat <- list()
ptm <- proc.time()
output.relat <- lapply(variables.r, simulation, model="model2",output="relatedness")  
# Run time (min)
cat(paste("Total runs:", length(variables.r), "; Simulation time:", round(((proc.time() - ptm)[3])/60, digits=2), "min"))

# Labelling model replicates
repliname=vector()
for(i in 1:length(variables)){
  ip = as.numeric(variables[[i]])
  repliname[i] <- paste("N",ip[1],'r',ip[2],'tprob',ip[3], sep="")
}
names(output.relat) <- repliname


# Organizing model replicates; every 100 replicates of a given set of params: output.relat.replic[[params]][[replic]][[time.step]][[output]]
output.relat.replic <- list(output.relat[c(1:100)], output.relat[c(101:200)], output.relat[c(201:300)], output.relat[c(301:400)], output.relat[c(401:500)], output.relat[c(501:600)], output.relat[c(601:700)], output.relat[c(701:800)])
names(output.relat.replic) <- names(output.relat)[c(100,200,300,400,500,600,700,800)]

# Picking one replicate of each parameter combination
params
output.relat.samples <- output.relat[c(50,150,250,350,450,550,650,750)]
names(output.relat.samples) <- names(output.relat)[c(100,200,300,400,500,600,700,800)]

# Saving data from simulation and all variables
save.image(file=paste(getwd(), "/data/3_relatedness_samples.RData", sep=""))










# 3. Relatedness Ratio  --------------------------------------------------------------


## Plotting 1 model run across the representative parameter space

# 3.1. Plot of Ratio and log-ratio of the average relatedness within foraging groups and the average relatedness expected by chance for each time of the simumation, across parameter space

# Permutation test for relatedness within foraging groups: Comparing the observed relatedness to the expected mean relatedness from the simulation data. First, we pick at random the same number of individuals in the same number of foraging groups, and calculate their mean relatedness. Do this randomisation test for 1000 times (iterations) for each of the time step of the simulation. Sencond, calculate the ratio of observed relatedness of the real foraging group to the mean random relatedness. Then plotting the log of this ratio give us information on how much the relatedness in the groups depart from what we expect by chance. 

# Plotting log relatedness ratio with random benchmark
par(mfrow=c(4,2), mar=c(4,4,2,0.5))
for(i in 1:nrow(params)){
  plot(1:n.reps, log(output.relat.samples[[i]][[5]]$all.relatedness$relat.ratio),las=1, ylab="", xlab="", main=paste("N =", params[i,1], ", r =", params[i,2], ', c =', params[i,3]))
  if(i %in% c(1,3,5,7)){ title(ylab="Log Mean Member Relatedness/Random")}
  if(i %in% 7:8 ){ title(xlab="Time step")}    
}

# Plotting relatedness ratio with random benchmark
par(mfrow=c(4,2), mar=c(4,4,2,0.5))
for(i in 1:nrow(params)){
  plot(1:n.reps, output.relat.samples[[i]][[5]]$all.relatedness$relat.ratio,las=1, ylab="", xlab="", main=paste("N =", params[i,1], ", r =", params[i,2], ', c =', params[i,3]))
  if(i %in% c(1,3,5,7)) { title(ylab="Mean Member Relatedness/Random")}
  if(i %in% 7:8){ title(xlab="Time step")}    
}





# 3.2 Plot of ratio and log-ratio of the average relatedness within foraging groups and all non-members individuals for each time of the simulation

# Plotting log relatedness ratio
par(mfrow=c(4,2), mar=c(4,4,2,0.5))
for(i in 1:nrow(params)){
  plot(1:n.reps, log(output.relat.samples[[i]][[5]]$all.relatedness$mean.coop.relatedness/output.relat.samples[[i]][[5]]$all.relatedness$mean.noncoop.relatedness),las=1, ylab="", xlab="", main=paste("N =", params[i,1], ", r =", params[i,2], ', c =', params[i,3]))
  if(i %in% c(1,3,5,7)){ title(ylab="Log Mean Relatedness Member/Non-member")}
  if(i %in% 7:8){ title(xlab="Time step")}    
}

# Plotting relatedness ratio
par(mfrow=c(4,2), mar=c(4,4,2,0.5))
for(i in 1:nrow(params)){
  plot(1:n.reps, output.relat.samples[[i]][[5]]$all.relatedness$mean.coop.relatedness/output.relat.samples[[i]][[5]]$all.relatedness$mean.noncoop.relatedness,las=1, ylab="", xlab="", main=paste("N =", params[i,1], ", r =", params[i,2], ', c =', params[i,3]))
  if(i %in% c(1,3,5,7)){ title(ylab="Mean Relatedness Member/Non-member")}
  if(i %in% 7:8 ){ title(xlab="Time step")}    
}






# 3.3 Plot of relatedness among foraging group members, non-members and random

# Plotting mean relatedness among members, among non-members, and random
par(mfrow=c(4,2), mar=c(4,4,2,0.5))
for(i in 1:nrow(params)){
  lim=round(max(output.relat.samples[[i]][[5]]$all.relatedness[,1:3], na.rm=T), digits=2)
  plot(1:n.reps, output.relat.samples[[i]][[5]]$all.relatedness$mean.random.relatedness, col=rgb(red=1,green=0,blue=0, alpha=0.1),ylab="", xlab="", las=1,main=paste("N =", params[i,1], ", r =", params[i,2], ', c =', params[i,3]),ylim=c(0,lim))
  points(1:n.reps, output.relat.samples[[i]][[5]]$all.relatedness$mean.noncoop.relatedness, col=rgb(red=0,green=0,blue=1, alpha=0.05),ylim=c(0,lim))
  points(1:n.reps, output.relat.samples[[i]][[5]]$all.relatedness$mean.coop.relatedness,ylim=c(0,lim))
  
  if(i == 2) legend('topright', c("Random", "Non-memb", "Memb"), cex=0.9, col=c("red","blue", "black"), bty="n", lwd=c(3, 3)) 
  if(i %in% c(1,3,5,7)){ title(ylab="Mean Relatedness")}
  if(i %in% 7:8 ){ title(xlab="Time step")}    
}


# Plotting relatedness among foraging group members and 95%CI random
par(mfrow=c(4,2), mar=c(4,4,2,0.5))
for(i in 1:nrow(params)){
  lim=round(max(output.relat.samples[[i]][[5]]$all.relatedness[,c(1,4,5)], na.rm=T), digits=2)
  
  errbar(x=1:n.reps, y=output.relat.samples[[i]][[5]]$all.relatedness$mean.coop.relatedness, output.relat.samples[[i]][[5]]$all.relatedness$low2.5CIrandom, output.relat.samples[[i]][[5]]$all.relatedness$upper97.5CIrandom, cap=0.005, col="black", ylab="", xlab="", colwhisker=rgb(red=1,green=0,blue=0, alpha=0.05),main=paste("N =", params[i,1], ", r =", params[i,2], ', c =', params[i,3]), ylim=c(0,lim))
  
  if(i == 2) {legend('topright', c("Random", "Memb"), cex=0.9, col=c("red", "black"), bty="n", lwd=c(3, 3))} 
  if(i %in% c(1,3,5,7)){ title(ylab="Mean member relatedness (95%CI)")}
  if(i %in% 7:8 ){ title(xlab="Time step")}  
}




## Considering all 100 model replicates

# 3.4 Members, non-members and random relatedness

# Plot across representative parameter space
par(mfrow=c(4,2), mar=c(4,4,2,0.5))
for(j in 1:nrow(params)){
  # relatedness y-axis limits
  lim=numeric()
  for(z in 1:length(output.relat.replic[[j]])){
    aux=output.relat.replic[[j]][[z]]
    lim[z]=max(aux[[5]]$all.relatedness[,1:3], na.rm=T)
  }
  # empty plot
  plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(0, round(max(lim),digits=2)), las=1, main=paste("N=", params[j,1], ", r=", params[j,2], ', c=', params[j,3], ', replic=100',sep=''))
  # overlay members, non-members, random relatedness
  for(i in 1:length(output.relat.replic[[j]])){
    # output.relat.replic[[params]][[replic]][[time.step]][[output]]
    points(1:n.reps, output.relat.replic[[j]][[i]][[5]]$all.relatedness$mean.random.relatedness, col=rgb(red=1,green=0,blue=0,alpha=0.01))
    points(1:n.reps, output.relat.replic[[j]][[i]][[5]]$all.relatedness$mean.noncoop.relatedness, col=rgb(red=0,green=0,blue=1,alpha=0.02))
    points(1:n.reps, output.relat.replic[[j]][[i]][[5]]$all.relatedness$mean.coop.relatedness, col=rgb(red=0,green=0,blue=0,alpha=0.3))
  }
  if(j == 2) legend('topright', c("Random", "Non-members", "Members"), cex=0.9, col=c("red","blue", "black"), bty="n", lwd=c(3, 3)) 
  if(j %in% c(1,3,5,7)){ title(ylab="Mean Relatedness")}
  if(j %in% 7:8 ){ title(xlab="Time step")}   
}




# 3.5 Log relatedness ratio

# Plotting log relatedness ratio
par(mfrow=c(4,2), mar=c(4,4,2,0.5))
for(j in 1:nrow(params)){
  # relatedness y-axis limits
  ulim <- llim <- numeric()
  for(z in 1:length(output.relat.replic[[j]])){
    aux=output.relat.replic[[j]][[z]]
    lim=log(aux[[5]]$all.relatedness$relat.ratio)
    lim[!is.finite(lim)] <- NA
    ulim[z]=max(lim, na.rm=T, is.finite=T)
    llim[z]=min(lim, na.rm=T, is.finite=T)     
  }
  # empty plot
  plot(1, type="n", xlab="", ylab="", xlim=c(0, n.reps), ylim=c(round(min(llim,na.rm=T),digits=1), round(max(ulim,na.rm=T),digits=1)), las=1, main=paste("N=", params[j,1], ", r=", params[j,2], ', c=', params[j,3], ', replic=100',sep=''))
  # overlay relatedness ratios
  for(i in 1:length(output.relat.replic[[j]])){
    points(1:n.reps, log(output.relat.replic[[j]][[i]][[5]]$all.relatedness$relat.ratio), col=rgb(red=0,green=0,blue=0,alpha=0.2))
  }
  if(j %in% c(1,3,5,7)){ title(ylab="Log Mean Relatedness Members/Random")}
  if(j %in% 7:8 ){ title(xlab="Time step")}   
}








# 4. Relatedness distribution --------------------------------------------------------------


#  Relatedness distribution: plots of relatedness among foraging group members and non-members

# Preparing data
for(i in 1:nrow(params)){
  relat.data <- vector("list", length(time.steps))
  for(which.sample in 1:5){
    relat <- output.relat.samples[[i]][[which.sample]]$relatedness.network[lower.tri(output.relat.samples[[i]][[which.sample]]$relatedness.network)]
    step <- rep(time.steps[which.sample], times=length(output.relat.samples[[i]][[which.sample]]$relatedness.network[lower.tri(output.relat.samples[[i]][[which.sample]]$relatedness.network)]))
    group <- output.relat.samples[[i]][[which.sample]]$coop.NOW[lower.tri(output.relat.samples[[i]][[which.sample]]$coop.NOW)]
    relat.data[[which.sample]] <- cbind(relat, step, group)
  }
  aux.data <- aux.data2 <- do.call(rbind, relat.data)
  aux.data <- aux.data2 <- as.data.frame(aux.data)
  colnames(aux.data) <- colnames(aux.data2) <- c("relatedness", "time.step", "groups")
  
  aux.data$time.step <- as.factor(aux.data$time.step)
  aux.data$groups <- as.factor(aux.data$groups)
  box.data[[i]] <- aux.data
  
  # data prep to overlay  plots. for foraging group members, time steps will be 300, 500, 700, 900, 1100
  aux.data2[which(aux.data2$time.step==200 & aux.data2$groups==1), ]$time.step = 300
  aux.data2[which(aux.data2$time.step==400 & aux.data2$groups==1), ]$time.step = 500
  aux.data2[which(aux.data2$time.step==600 & aux.data2$groups==1), ]$time.step = 700
  aux.data2[which(aux.data2$time.step==800 & aux.data2$groups==1), ]$time.step = 900
  aux.data2[which(aux.data2$time.step==1000 & aux.data2$groups==1), ]$time.step = 1100
  aux.data2$time.step <- as.factor(aux.data2$time.step)
  aux.data2$groups <- as.factor(aux.data2$groups)
  box.data2[[i]] <- aux.data2
}

# Plotting relatedness: foraging group
pl<-list()
for(i in 1:nrow(params)){
  g1=box.data[[i]][which(box.data[[i]]$groups==1),]
  pl[[i]] <-
    ggplot(data=g1, aes(x=time.step,y=relatedness)) +
    geom_jitter(shape=1, size=1, position=position_jitter(0.5, 0.1), col="black") +
    geom_violin(fill=NA,colour="grey", size=1, scale="count") +
    stat_summary(fun.data="mean_sdl", fun.args=list(mult=1), geom="pointrange", color="red", size=1, shape=20)+ 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.text.x=element_text(size=14), 
          axis.text.y=element_text(size=14), 
          axis.title=element_text(size=14)) +
    scale_y_continuous(limits = c(-0.2,0.6)) +
    labs(x="Time step",y="Pairwise relatedness") +
    ggtitle(paste("N =", params[i,1], ", r =", params[i,2], ', c =', params[i,3]))
}
grid.arrange(pl[[1]], pl[[2]],pl[[3]],pl[[4]],
             pl[[5]], pl[[6]],pl[[7]],pl[[8]],
             pl[[9]], pl[[10]],pl[[11]],pl[[12]],nrow=3)


# Plotting relatedness: foraging group members and non-members
pl2<-list()
for(i in 1:nrow(params)){  
  o <- ggplot(box.data2[[i]], aes(x=time.step, y=relatedness, col=groups))
  pl2[[i]] <-
    o + geom_jitter(shape=1, size=0.1, position=position_jitter(1)) +
    stat_summary(fun.data="mean_sdl", fun.args=list(mult=1), geom="pointrange", aes(color=paste(groups), group=groups), size=1, shape=20, show.legend=NA) +
    theme(legend.position=c(0.10,0.70), 
          legend.title=element_blank(),
          legend.background = element_rect(colour = "transparent", fill = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.text.x=element_text(size=14), 
          axis.text.y=element_text(size=14), 
          axis.title=element_text(size=14)) +
    scale_color_manual(labels = c("non", "memb"), values = c("grey", "red")) +
    scale_x_discrete("Time Step", breaks=c('200', '400','600','800', '1000')) +   
    scale_y_continuous("Relatedness") +
    ggtitle(paste("N =", params[i,1], ", r =", params[i,2], ', c =', params[i,3]))
}
grid.arrange(pl2[[1]], pl2[[2]],pl2[[3]],pl2[[4]],
             pl2[[5]], pl2[[6]],pl2[[7]],pl2[[8]],
             pl2[[9]], pl2[[10]],pl2[[11]],pl2[[12]],nrow=3)










# 5. Network Plots ----------------------------------------------------------------

# Refreshing the representative parameter space (Add edge filter for social network)
N <- c(40,120)          # Population/network size
r <- c(5,15,25,35)      # Resource patch size
tprob <- c(0.2)         # Initial network connectance
n.reps <- 1000          # Time steps
replic <- 100           # Number of model replicates
time.steps=seq(from=(n.reps/5), to=n.reps, by=n.reps/5)# Output samples

params <- expand.grid(N,r,tprob)
colnames(params) <- c("N","r","tprob")

params <- cbind(params, rep(c(0.25, 0.45), nrow(params)/2)); colnames(params) <- c('N','r','tprob','filter')

# Refreshing the picked up simulating samples across the parameter space
output.relat.samples <- output.relat[c(100,200,300,400,500,600,700,800)]
names(output.relat.samples) <- names(output.relat)[c(100,200,300,400,500,600,700,800)]



# 5.1 Pedigree, relatedness and social network, identifying foraging group members

# To select a given set of parameters: which.net='X'; params[which.net,]

for(which.net in 1:nrow(params)){
  
  par(mfcol=c(3,5), mar=c(0.5,4,1,0.5))
  
  for (which.sample in 1:5) {   
    
    # Plot pedigree network
    # Alive individuals (members=large red nodes, non-members=small black nodes) connected by pedigree
    graph2 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.link,mode="undirected",diag=FALSE,weighted=TRUE)
    plot(graph2,
         edge.width=(exp(E(graph2)$weight)-1),
         vertex.size=(as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1)*3,
         vertex.label="",
         mark.groups=output.relat.samples[[which.net]][[which.sample]]$groups, 
         vertex.color=c("black","red")[as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1],
         main=paste('Time step = ', time.steps[which.sample]))
    if(which.sample == 1){ title(ylab="Pedigree", font=2)}  
    
    #  plot relatedness network
    graph3 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.network,mode="undirected",diag=FALSE,weighted=TRUE)
    plot(graph3,
         edge.width=10*(exp(E(graph3)$weight)-1)^2,
         vertex.size=(as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1)*3,
         vertex.label="",
         vertex.color=c("black","red")[as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1],
         mark.groups=output.relat.samples[[which.net]][[which.sample]]$groups) 
    if(which.sample == 1){ title(ylab="Relatedness", font=2)}
    
    # Plot social network
    # Alive individuals connected by proportion of times seen in foraging groups during the simulation time (filtered links)
    # Creating social network links as the proportion of times individuals cooperated during the model run
    
    # Creating social network links as the proportion of times ALIVE individuals cooperated during the model run
    N.inds.alive <- sum(output.relat.samples[[which.net]][[which.sample]]$ids$Alive==1)
    coop.total.tmp <- matrix(0,N.inds.alive,N.inds.alive)
    for (i in 1:(N.inds.alive-1)) {
      for (j in (i+1):N.inds.alive) {
        if (sum(output.relat.samples[[which.net]][[which.sample]]$alive[,i]==1 & output.relat.samples[[which.net]][[which.sample]]$alive[,j]==1)>0) {
          coop.total.tmp[i,j] <- output.relat.samples[[which.net]][[which.sample]]$coop.total[i,j]/sum(output.relat.samples[[which.net]][[which.sample]]$alive[,i]==1 & output.relat.samples[[which.net]][[which.sample]]$alive[,j]==1)
          coop.total.tmp[j,i] <- coop.total.tmp[i,j]
        }
      }
    }
    # filtering network edges for better visualization
    edge.filter = params[which.net,4]
    coop.total.tmp2 <- coop.total.tmp
    coop.total.tmp2[coop.total.tmp2<edge.filter] <- 0
    
    graph4 <- graph.adjacency(coop.total.tmp2,mode="undirected",diag=FALSE,weighted=TRUE)
    plot(graph4,
         edge.width=(exp(E(graph4)$weight)-1),
         vertex.size=(as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1)*3,
         vertex.label="",
         mark.groups=output.relat.samples[[which.net]][[which.sample]]$groups, 
         vertex.color=c("black","red")[as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1])
    
    if(which.sample == 1){ title(ylab="Social", font=2)}
    
  } 
  
}







# 5.2. Full pedigree, full relatedness, full social network


for(which.net in 1:nrow(params)){
  
  par(mfcol=c(3,5),mar=c(0.5,4,1,0.5))
  
  for (which.sample in 1:5) {
    
    # Full pedigree network
    # All individuals (dead=black, alive=green) connected by pedigree during the simulation time
    graph5 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.link.total,mode="undirected",diag=FALSE,weighted=TRUE)
    plot(graph5,
         edge.width=(exp(E(graph5)$weight)-1),
         vertex.size=5,
         vertex.label="",
         vertex.color=c("black","green")[output.relat.samples[[which.net]][[which.sample]]$ids$Alive+1],
         main=paste('Time step = ', time.steps[which.sample]))
    if(which.sample == 1){ title(ylab="Pedigree", font=2)}    
    
    # Full relatedness network
    # All individuals (dead=black, alive=green) connected by total relatedness during the simulation time
    graph6 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.network.total,mode="undirected",diag=FALSE,weighted=TRUE)
    plot(graph6,
         edge.width=10*(exp(E(graph6)$weight)-1)^2,
         vertex.size=5,
         vertex.label="",
         vertex.color=c("black","green")[output.relat.samples[[which.net]][[which.sample]]$ids$Alive+1])
    if(which.sample == 1){ title(ylab="Relatedness", font=2)}    
    
    
    # Full social network
    # All individuals (dead=black, alive=green) connected by proportion of times seen in foraging groups during the simulation time (filtered links <0.05)
    # Creating social network links as the proportion of times individuals cooperated during the model run
    N.inds <- nrow(output.relat.samples[[which.net]][[which.sample]]$ids)
    coop.total.tmp <- matrix(0,N.inds,N.inds)
    for (i in 1:(N.inds-1)) {
      for (j in (i+1):N.inds) {
        if (sum(output.relat.samples[[which.net]][[which.sample]]$alive[,i]==1 & output.relat.samples[[which.net]][[which.sample]]$alive[,j]==1)>0) {
          coop.total.tmp[i,j] <- output.relat.samples[[which.net]][[which.sample]]$coop.total[i,j]/sum(output.relat.samples[[which.net]][[which.sample]]$alive[,i]==1 & output.relat.samples[[which.net]][[which.sample]]$alive[,j]==1)
          coop.total.tmp[j,i] <- coop.total.tmp[i,j]
        }
      }
    }
    # filtering network edges for better visualization
    edge.filter = params[which.net,4]
    coop.total.tmp2 <- coop.total.tmp
    coop.total.tmp2[coop.total.tmp2<edge.filter] <- 0
    # Making graph and plotting social network
    graph7 <- graph.adjacency(coop.total.tmp2,mode="undirected",diag=FALSE,weighted=TRUE)
    plot(graph7,
         edge.width=(exp(E(graph7)$weight)-1),
         vertex.size=5,
         vertex.label="",
         vertex.color=c("black","green")[output.relat.samples[[which.net]][[which.sample]]$ids$Alive+1])
    if(which.sample == 1){ title(ylab="Social", font=2)}    
  }
  
}






# 5.3. Social network with clustering on relatedness
# colour nodes by family groups: run the clustering algorithm on the relatedness network, but plot the exclusivity network

for(which.net in 1:nrow(params)){
  
  par(mfcol=c(1,5), mar=c(0.5,4,1,0.5))
  
  for (which.sample in 1:5) {
    
    # relatedness network with clustering
    graph8 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.network,mode="undirected",diag=FALSE,weighted=TRUE)
    modul <- cluster_walktrap(graph8, weights = E(graph8)$weight, steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE)
    
    # Creating social network links as the proportion of times ALIVE individuals cooperated during the model run
    N.inds.alive <- sum(output.relat.samples[[which.net]][[which.sample]]$ids$Alive==1)
    coop.total.tmp <- matrix(0,N.inds.alive,N.inds.alive)
    for (i in 1:(N.inds.alive-1)) {
      for (j in (i+1):N.inds.alive) {
        if (sum(output.relat.samples[[which.net]][[which.sample]]$alive[,i]==1 & output.relat.samples[[which.net]][[which.sample]]$alive[,j]==1)>0) {
          coop.total.tmp[i,j] <- output.relat.samples[[which.net]][[which.sample]]$coop.total[i,j]/sum(output.relat.samples[[which.net]][[which.sample]]$alive[,i]==1 & output.relat.samples[[which.net]][[which.sample]]$alive[,j]==1)
          coop.total.tmp[j,i] <- coop.total.tmp[i,j]
        }
      }
    }
    # filtering network edges for better visualization
    edge.filter = params[which.net,4]
    coop.total.tmp2 <- coop.total.tmp
    coop.total.tmp2[coop.total.tmp2<edge.filter] <- 0
    
    # plot social network at the given time step with clustering from relatedness
    graph9 <- graph.adjacency(coop.total.tmp2,mode="undirected",diag=FALSE,weighted=TRUE)
    
    V(graph9)$color <- modul$membership + 1
    
    plot(graph9,
         edge.width=(exp(E(graph9)$weight)),
         vertex.size=(as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1)*3,
         vertex.label="",
         mark.groups=output.relat.samples[[which.net]][[which.sample]]$groups,
         main=paste('Time step = ', time.steps[which.sample]))
  }
  
}








# 6. Generations ----------------------------------------------------------

# Calculating the average number of generations. The diameter (i.e. longest path) of the Full Pedigree Network at time step=1000 is the maximum number of generations in given model run. This value is then averaged across all 100 model replicates for a given parameter combination

maxgen <- ngroup <- mgroup <- rep( list(numeric()), nrow(params) )
meangen <- numeric()
for(j in 1:nrow(params)){
  for(i in 1:length(output.relat.replic[[j]])){
    # output.relat.replic[[params]][[replic]][[time.step]][[output]]
    graph5 <- graph.adjacency(output.relat.replic[[j]][[i]][[5]]$relatedness.link.total,mode="undirected",diag=FALSE,weighted=TRUE)
    # calculate max number of generations
    maxgen[[j]][i] <- diameter(graph5, directed = F, unconnected = TRUE, weights = NULL)
    # calculate number of groups
    ngroup[[j]][i] <- length(output.relat.replic[[j]][[i]][[5]]$groups)
    # calculate mean group size
    mgroup[[j]][i] <- mean(unlist(lapply(output.relat.replic[[j]][[i]][[5]]$groups, length)))
  }
  meangen[j] <- mean(maxgen[[j]])
}

# mean number of generations across parameter space
data.frame(params, 'mean.generation'=meangen)

# Overall average number of generations
mean(meangen)


# plot generation, mean group size, number of groups
par(mfrow=c(2,4))
dataplot <- list()
for(j in 1:nrow(params)) { 
  dataplot[[j]] <- cbind(unlist(maxgen[[j]]), unlist(ngroup[[j]]), unlist(mgroup[[j]]))
  plot(jitter(dataplot[[j]][,1]), (dataplot[[j]][,3]), xlab="generation", ylab="groups", ylim=c(0,40), xlim=c(1,9))
  points(jitter(dataplot[[j]][,1]), (dataplot[[j]][,2]), col='red')
  if(j==1) {legend(2, 35, legend=c("mean size", "number"), fill=c("black", "red"), cex=0.8)}
  mtext(paste("N=",params[j,1],', r=',params[j,2], ', Tprob=', params[j,3], sep=""), side=3)
}
