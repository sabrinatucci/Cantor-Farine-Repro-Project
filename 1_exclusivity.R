# Exclusivity of foraging group

# Summary -----------------------------------------------------------------
#
# 1. Initial setup
# 2. Model run
# 3. Results: surface plot
# 4. Results: network plot
# 5. Saving working space



# 1. Initial setup --------------------------------------------------------

# Load packages and functions
source("setup.R")

# Define parameter space
N <- c(20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200)  # Population/network size
r <- c(5,7,9,11,13,15,17,21,25,29,31,37,41,45,51)                 # Resource patch size
tprob <- c(0.2,0.5,0.8)                                           # Initial network connectance
replic <- 100                                                     # Number of model replicates
payoff.type <- 'size-based'                                       # Type of group resource share

variables <- expand.grid(N,r,tprob,payoff.type)
colnames(variables) <- c("N","r","tprob","type")
variables <- variables[which(variables$r<=(variables$N/2)),]      # only when resource.size <= half of population
variables <- variables[rep(seq_len(nrow(variables)), each=replic),]
variables <- split(variables, seq(nrow(variables)))



# 2. Model run ------------------------------------------------------------

# Run model1 across the parameter space, and outputs a matrix with 4 columns ("N","r","tprob","Exclusivity") for the surface plot
ptm <- proc.time()
output <- do.call('rbind',lapply(variables, simulation, model="model1", output="exclusivity"))
colnames(output) <- c("N","r","tprob","Exclusivity")
cat(paste("simulation time:", round(((proc.time() - ptm)[3])/60, digits=2), "min"))


# Alternatively: load simulation data
#load(paste(getwd(), "/data/1_exclusivity_surface.RData", sep=""))
# only tprob=0.5
#load(paste(getwd(), "/data/1_exclusivity_surface_tprob05.RData", sep=""))




# 3. Results: surface plot --------------------------------------------------

# Preparing the data for surface plot
output <- output[which(!is.na(output[,4])),]
output <- as.data.frame(output)
output$ID <- 1:nrow(output)

# Bins for the plot
r.bins <- seq(min(output$r),max(output$r),1)
N.bins <- seq(min(output$N),max(output$N),1)

par(mfrow=c(1,3), mar=c(1,1,1,1))
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
    res <- persp(r.bins, N.bins, matrix(NA,nrow=length(r.bins),ncol=length(N.bins)), theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors[facetcol], ylab="",xlab="",zlab="",zlim=c(zmin,zmax),xlim=range(r.bins),ylim=range(N.bins),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75, main=paste("Connectance =", tprob[s]))
    
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
    res <- persp(r.bins, N.bins, z, theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors2[facetcol], ylab="Population size",xlab="Resource patch size",zlab="Exclusivity",zlim=c(zmin,zmax),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75)
    input3 <- input2[which((input2$Exclusivity - input2$pred) >= 0),]
    for (i in 1:nrow(input3)) {
      lines(trans3d(c(input3$r[i],input3$r[i]),c(input3$N[i],input3$N[i]),c(input3$Exclusivity[i],input3$pred[i]),res), pch=20,col="grey")
    }
    points(trans3d(input3$r,input3$N,input3$Exclusivity,res), pch=20)

}


# 4. Results: networks plot -----------------------------------------------

# plotting a typical social network created from model 1

# Typical model run
network1 = model1(N=40, resource.size=5, n.reps=100, tprob=0.2, type="size-based", output="model")
  
# Creating social network links as the proportion of times individuals cooperated during the model run
final.adj.matrix = network1$coop.total
n.reps = 100

# then filtering network edges for better visualization
edge.filter = 0.15
coop.total.tmp <- final.adj.matrix/n.reps
coop.total.tmp[coop.total.tmp<edge.filter] <- 0

# Making graph and plotting
graph1 <- graph.adjacency(coop.total.tmp,mode="undirected",diag=FALSE,weighted=TRUE)
par(mfrow=c(1,1),mar=c(1,1,1,1))
plot(graph1,edge.width=(exp(E(graph1)$weight))*1,vertex.size=10,vertex.color='grey',vertex.label="", layout=layout_nicely(graph1))



# Plotting network examples across parameter space
par(mfrow=c(2,2), mar=c(1,1,1,1))
params <- rbind(c(80,10,0.4), 
                c(80,40,0.4), 
                c(180,10,0.4), 
                c(180,40,0.4))
colnames(params) <- c('N','r','filter')

for(i in 1:nrow(params)){
    net = model1(N=params[i,1], resource.size=params[i,2], n.reps=100, tprob=0.2, type="size-based", output="model")  
    final.adj.matrix = net$coop.total
    n.reps = 100
    edge.filter = params[i,3]
    coop.total.tmp <- final.adj.matrix/n.reps
    coop.total.tmp[coop.total.tmp<edge.filter] <- 0
    g <- graph.adjacency(coop.total.tmp,mode="undirected",diag=FALSE,weighted=TRUE)
    plot(g,
         edge.width=(exp(E(g)$weight))*1.2,
         vertex.color='grey',
         edge.color="black",
         vertex.label="",
         vertex.size=(as.numeric(rowSums(net$coop.NOW)>=1)+1)*3,
#         vertex.color=c("white","red")[as.numeric(rowSums(net$coop.NOW)>=1)+1],
         mark.groups=net$groups,
         main=paste("N =", params[i,1], ", r =", params[i,2], ', filter =', params[i,3]))    
}





# 5. Saving working space  ---------------------------------------------------------------

# saving data from simulation and all variables
#save.image(file=paste(getwd(), "/data/1_exclusivity_surface.RData", sep="")) #will save anything in environment 