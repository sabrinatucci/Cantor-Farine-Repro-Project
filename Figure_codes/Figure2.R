# Figure 2: exclusivity surface plot and sample social networks across the representative parameter space

# Figure 2 -------------------------------------------------------------

# Load data and functions
load(paste(getwd(), "/data/1_exclusivity_surface.RData", sep=""))
source("setup.R")

layout(matrix(c(0,2,0,
                   3,1,4,
                   0,5,0),3,3,byrow = TRUE), c(1,2), c(1,2), F)


# Exclusivity Surface plot for connectance = 0.2

output <- output[which(!is.na(output[,4])),]
output <- as.data.frame(output)
output$ID <- 1:nrow(output)

r.bins <- seq(min(output$r),max(output$r),1)
N.bins <- seq(min(output$N),max(output$N),1)

input <- output[which(output$tprob==0.2 & output$r<=(output$N/2)),]

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
res <- persp(r.bins, N.bins, matrix(NA,nrow=length(r.bins),ncol=length(N.bins)), theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors[facetcol], ylab="",xlab="",zlab="",zlim=c(zmin,zmax),xlim=range(r.bins),ylim=range(N.bins),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75)

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




# Plotting sample network across the representative parameter space
par(mfrow=c(2,2), mar=c(1,1,1,1))
params <- rbind(c(60,10,0.3), 
                c(60,30,0.3), 
                c(160,10,0.4), 
                c(160,30,0.4))
net <- list()
for(i in 1:nrow(params)){
  net[[i]] = model1(N=params[i,1], resource.size=params[i,2], n.reps=100, tprob=0.2, output="model")  
  final.adj.matrix = net[[i]]$coop.total
  n.reps = 100
  coop.total.tmp <- final.adj.matrix/n.reps
  
  #gaux <- graph.adjacency(coop.total.tmp,mode="undirected",diag=FALSE,weighted=TRUE)
  #coords <- layout_with_fr(gaux,weights=E(gaux)$weight)
  
  edge.filter = params[i,3]
  coop.total.tmp[coop.total.tmp<edge.filter] <- 0
  
  # Define coordinates for plotting group members nodes closer
  social.network <- coop.total.tmp
  N <- nrow(social.network)
  group.net <- matrix(0,N,N)
  groups <- net[[i]]$groups
  for (j in 1:length(groups)) {
    group.net[groups[[j]],groups[[j]]] <- sample(c(0,1),length(groups[[j]])^2,prob=c(0.9,0.1),replace=TRUE)
  }
  composite.network <- 100000*group.net + 1000*social.network
  graph.composite.network <- graph.adjacency(composite.network,mode="undirected",diag=FALSE,weighted=TRUE)
  coords <- layout_with_fr(graph.composite.network,weights=E(graph.composite.network)$weight)
  
  g <- graph.adjacency(coop.total.tmp,mode="undirected",diag=FALSE,weighted=TRUE)
  

  plot(g,
       layout=coords,
       edge.width=(exp(E(g)$weight)),
       edge.color="grey",
       vertex.size=(as.numeric(rowSums(net[[i]]$coop.NOW)>=1)+1)*3,
       vertex.color=c("black","red")[as.numeric(rowSums(net[[i]]$coop.NOW)>=1)+1],
       vertex.label="",
       vertex.shape="fcircle", 
       vertex.frame.color="black",
       vertex.frame.width=0.5
       #       mark.groups=net[[i]]$groups
       )
}
