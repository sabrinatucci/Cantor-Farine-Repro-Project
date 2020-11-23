# Figure S3: networks across parameter space --------------------------------------------------------

# Refreshing the representative parameter space (Add edge filter for social network)
N <- c(40,120)          # Population/network size
r <- c(5,15,25,35)      # Resource patch size
tprob <- c(0.2)         # Initial network connectance
n.reps <- 1000          # Time steps
replic <- 100           # Number of model replicates
time.steps=seq(from=(n.reps/5), to=n.reps, by=n.reps/5)# Output samples

params <- expand.grid(N,r,tprob)
colnames(params) <- c("N","r","tprob")

# Load packages and data
source("setup.R")
load(paste(getwd(), "/data/3_relatedness_samples.RData", sep=""))

params <- cbind(params, rep(c(0.3, 0.6), nrow(params)/2)); colnames(params) <- c('N','r','tprob','filter')

# Picking up simulating samples across the parameter space
params
output.relat.samples <- output.relat[c(50,150,250,350,450,550,650,750)]
names(output.relat.samples) <- names(output.relat)[c(100,200,300,400,500,600,700,800)]

#Choose which parameter combination to plot
params[c(3,7,4,8),]


# Plotting networks
for(which.net in c(3,7,4,8)) {
  
  par(mfcol=c(6,5), mar = c(0.1, 3, 1, 0.1), oma = c(0.1, 0.1, 0.1, 0.1), tcl = -0.25, mgp = c(2, 0.6, 0))
     
    for (which.sample in 1:5) {   
      
      ### (a) Pedigree network
      
      ## (a1) Alive individuals (foraging group members=large red nodes, non-members=small black nodes) connected by pedigree
      graph2 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.link,mode="undirected",diag=FALSE,weighted=TRUE)
      
      # Define coordinates for plotting group members nodes closer
      relatedness.link <- output.relat.samples[[which.net]][[which.sample]]$relatedness.link
      N <- nrow(relatedness.link)
      group.net <- matrix(0,N,N)
      groups <- output.relat.samples[[which.net]][[which.sample]]$groups
      for (i in 1:length(groups)) {
        group.net[groups[[i]],groups[[i]]] <- sample(c(0,1),length(groups[[i]])^2,prob=c(0.9,0.1),replace=TRUE)
      }
      composite.network <- 1000*group.net + 1000*relatedness.link
      graph.composite.network <- graph.adjacency(composite.network,mode="undirected",diag=FALSE,weighted=TRUE)
      coords <- layout_with_fr(graph.composite.network,weights=E(graph.composite.network)$weight)
      
      plot(graph2,
           edge.width=(exp(E(graph2)$weight)-1),
           vertex.size=(as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1)*3,
           vertex.label="",
           vertex.shape="fcircle", 
           vertex.frame.color="black",
           vertex.frame.width=0.5,
           mark.groups=output.relat.samples[[which.net]][[which.sample]]$groups, 
           vertex.color=c("black","red")[as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1],
           main=paste('Time step = ', time.steps[which.sample]),
           layout=coords)
      if(which.sample == 1){ title(ylab="(a) Pedigree", font=2)}  
      
      
      ## (a2) Full pedigree network: All individuals (dead=black, alive=green) connected by all pedigree during the simulation time
      # Identifying group members at each time step and update their indices, considering IDs of all agents (alive+dead) 
      alive.index = which(output.relat.samples[[which.net]][[which.sample]]$ids$Alive==1)
      coop.alive.index = alive.index[unlist(output.relat.samples[[which.net]][[which.sample]]$groups)]
      auxcoop=rep(0,nrow(output.relat.samples[[which.net]][[which.sample]]$ids))
      auxcoop[coop.alive.index]=1
      
      graph5 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.link.total,mode="undirected",diag=FALSE,weighted=TRUE)
      
      # Define coordinates for plotting group members nodes closer
      relatedness.link <- output.relat.samples[[which.net]][[which.sample]]$relatedness.link.total
      N <- nrow(relatedness.link)
      group.net <- matrix(0,N,N)
      #groups <- output.relat.samples[[which.net]][[which.sample]]$groups
      groups <- list(coop.alive.index)
      for (i in 1:length(groups)) {
        group.net[groups[[i]],groups[[i]]] <- sample(c(0,1),length(groups[[i]])^2,prob=c(0.9,0.1),replace=TRUE)
      }
      composite.network <- 10000*group.net + 10000*relatedness.link
      graph.composite.network <- graph.adjacency(composite.network,mode="undirected",diag=FALSE,weighted=TRUE)
      coords <- layout_with_fr(graph.composite.network,weights=E(graph.composite.network)$weight)
      
      plot(graph5,
           edge.width=(exp(E(graph5)$weight)-1),
           #vertex.size=5,
           vertex.size=(auxcoop+1)*3, # larger nodes represent alive group members
           vertex.shape="fcircle", 
           vertex.frame.color="black",
           vertex.frame.width=0.5,
           vertex.label="",
           vertex.color=c("black","green")[output.relat.samples[[which.net]][[which.sample]]$ids$Alive+1],
           #main=paste('Time step = ', time.steps[which.sample]),
           layout=coords)
      #if(which.sample == 1){ title(ylab="Pedigree", font=2)}    
      
      
      
      ### (b) Relatedness network
      
      ## (b1) Alive individuals (foraging group members=large red nodes, non-members=small black nodes) connected by relatedness
      graph3 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.network,mode="undirected",diag=FALSE,weighted=TRUE)
      
      # Define coordinates for plotting group members nodes closer
      relatedness.network <- output.relat.samples[[which.net]][[which.sample]]$relatedness.network
      N <- nrow(relatedness.network)
      group.net <- matrix(0,N,N)
      groups <- output.relat.samples[[which.net]][[which.sample]]$groups
      for (i in 1:length(groups)) {
        group.net[groups[[i]],groups[[i]]] <- sample(c(0,1),length(groups[[i]])^2,prob=c(0.9,0.1),replace=TRUE)
      }
      composite.network <- 1000*group.net + 1000*relatedness.network
      graph.composite.network <- graph.adjacency(composite.network,mode="undirected",diag=FALSE,weighted=TRUE)
      coords <- layout_with_fr(graph.composite.network,weights=E(graph.composite.network)$weight)
      
      plot(graph3,
           edge.width=3*(exp(E(graph3)$weight)-1)^2,
           vertex.size=(as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1)*3,
           vertex.label="",
           vertex.shape="fcircle", 
           vertex.frame.color="black",
           vertex.frame.width=0.5,
           vertex.color=c("black","red")[as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1],
           mark.groups=output.relat.samples[[which.net]][[which.sample]]$groups,
           layout=coords) 
      if(which.sample == 1){ title(ylab="(b) Relatedness", font=2)}
      
      
      ## (b2) Full relatedness network: All individuals (dead=black, alive=green) connected by total relatedness during the simulation time
      graph6 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.network.total,mode="undirected",diag=FALSE,weighted=TRUE)
      
      # Define coordinates for plotting group members nodes closer
      relatedness.network <- output.relat.samples[[which.net]][[which.sample]]$relatedness.network.total
      N <- nrow(relatedness.network)
      group.net <- matrix(0,N,N)
      #groups <- output.relat.samples[[which.net]][[which.sample]]$groups
      groups <- list(coop.alive.index)
      for (i in 1:length(groups)) {
        group.net[groups[[i]],groups[[i]]] <- sample(c(0,1),length(groups[[i]])^2,prob=c(0.9,0.1),replace=TRUE)
      }
      composite.network <- 10000*group.net + 10000*relatedness.network
      graph.composite.network <- graph.adjacency(composite.network,mode="undirected",diag=FALSE,weighted=TRUE)
      coords <- layout_with_fr(graph.composite.network,weights=E(graph.composite.network)$weight)
      
      plot(graph6,
           edge.width=3*(exp(E(graph6)$weight)-1)^2,
           #vertex.size=5,
           vertex.size=(auxcoop+1)*3, # larger nodes represent alive group members
           vertex.shape="fcircle", 
           vertex.frame.color="black",
           vertex.frame.width=0.5,
           vertex.label="",
           vertex.color=c("black","green")[output.relat.samples[[which.net]][[which.sample]]$ids$Alive+1],
           layout=coords)
      #if(which.sample == 1){ title(ylab="Relatedness", font=2)}   
      
      
      
      
      ### (c) Social networks
      
      ## (c1) Alive individuals (foraging group members=large red nodes, non-members=small black nodes) connected by gambit of the group using relatedness network to define clusters
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
      
      # Define coordinates for plotting group members nodes closer
      social.network <- coop.total.tmp
      N <- nrow(social.network)
      group.net <- matrix(0,N,N)
      groups <- output.relat.samples[[which.net]][[which.sample]]$groups
      for (i in 1:length(groups)) {
        group.net[groups[[i]],groups[[i]]] <- sample(c(0,1),length(groups[[i]])^2,prob=c(0.9,0.1),replace=TRUE)
      }
      composite.network <- 10000*group.net + 10000*social.network
      graph.composite.network <- graph.adjacency(composite.network,mode="undirected",diag=FALSE,weighted=TRUE)
      coords <- layout_with_fr(graph.composite.network,weights=E(graph.composite.network)$weight)
      
      # filtering network edges for better visualization
      edge.filter = params[which.net,4]
      coop.total.tmp2 <- coop.total.tmp
      coop.total.tmp2[coop.total.tmp2<edge.filter] <- 0
      
      # plot social network at the given time step with clustering (colors) from relatedness
      graph9 <- graph.adjacency(coop.total.tmp2,mode="undirected",diag=FALSE,weighted=TRUE)
      
      com=modul$membership
      comcolor <- colorRampPalette(c("white", "yellow", "red"), bias=10)( length(groups(modul)) )
      V(graph9)$color <- comcolor[com]
      
      plot(graph9,
           edge.width=(exp(E(graph9)$weight)),
           vertex.shape="fcircle", 
           vertex.frame.color="black",
           vertex.frame.width=0.5,
           vertex.size=(as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1)*5,
           vertex.label="",
           #mark.groups=output.relat.samples[[which.net]][[which.sample]]$groups,
           layout=coords)
      
      if(which.sample == 1){ title(ylab="(c) Social", font=2)}    
      
      
      ## (c2) Full social network: All individuals (dead=black, alive=green) connected by proportion of times seen in foraging groups during the simulation time
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
      
      # Define coordinates for plotting group members nodes closer
      social.network <- coop.total.tmp
      N <- nrow(social.network)
      group.net <- matrix(0,N,N)
      #groups <- output.relat.samples[[which.net]][[which.sample]]$groups
      groups <- list(coop.alive.index)
      for (i in 1:length(groups)) {
        group.net[groups[[i]],groups[[i]]] <- sample(c(0,1),length(groups[[i]])^2,prob=c(0.9,0.1),replace=TRUE)
      }
      composite.network <- 10000*group.net + 10000*social.network
      graph.composite.network <- graph.adjacency(composite.network,mode="undirected",diag=FALSE,weighted=TRUE)
      coords <- layout_with_fr(graph.composite.network,weights=E(graph.composite.network)$weight)
      
      # filtering network edges for better visualization
      edge.filter = params[which.net,4]
      coop.total.tmp2 <- coop.total.tmp
      coop.total.tmp2[coop.total.tmp2<edge.filter] <- 0
      
      # Making graph and plotting social network
      graph7 <- graph.adjacency(coop.total.tmp2,mode="undirected",diag=FALSE,weighted=TRUE)
      plot(graph7,
           edge.width=(exp(E(graph7)$weight)-1),
           #vertex.size=5,
           vertex.size=(auxcoop+1)*3, # larger nodes represent alive group members
           vertex.shape="fcircle", 
           vertex.frame.color="black",
           vertex.frame.width=0.5,
           vertex.label="",
           vertex.color=c("black","green")[output.relat.samples[[which.net]][[which.sample]]$ids$Alive+1],
           layout=coords)
      #if(which.sample == 1){ title(ylab="Social", font=2)} 
      
    } 
    
  }