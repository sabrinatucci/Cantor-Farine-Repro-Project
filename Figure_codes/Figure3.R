# Figure 3: networks --------------------------------------------------------

# Load packages and data
source("setup.R")
load(paste(getwd(), "/data/3_relatedness_samples.RData", sep=""))

# Refreshing the representative parameter space (Add edge filter for social network)
N <- c(40,120)          # Population/network size
r <- c(5,15,25,35)      # Resource patch size
tprob <- c(0.2)         # Initial network connectance
n.reps <- 1000          # Time steps
replic <- 100           # Number of model replicates
time.steps=seq(from=(n.reps/5), to=n.reps, by=n.reps/5)# Output samples

params <- expand.grid(N,r,tprob)
colnames(params) <- c("N","r","tprob")
params <- cbind(params, rep(c(0.3, 0.3), nrow(params)/2)); colnames(params) <- c('N','r','tprob','filter')

# Picking up simulating samples across the parameter space
params
output.relat.samples <- output.relat[c(51,151,251,351,451,551,651,751)]
names(output.relat.samples) <- names(output.relat)[c(c(51,151,251,351,451,551,651,751))]

#Choose which parameter combination to plot, using which.net=X, where X=nrow in params

which.net=3; params[which.net,]

# Plot networks ----
 par(mfcol=c(3,5), mar = c(0.1, 3, 1, 0.1), oma = c(0.1, 0.1, 0.1, 0.1), tcl = -0.25, mgp = c(2, 0.6, 0))

for (which.sample in 1:5) {   
    
    ## Plot pedigree network
    # Alive individuals (foraging group members=large red nodes, non-members=small black nodes) connected by pedigree
    graph2 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.link,mode="undirected",diag=FALSE,weighted=TRUE)
    
    # Define coordinates for plotting group members nodes closer
    relatedness.link <- output.relat.samples[[which.net]][[which.sample]]$relatedness.link
    N <- nrow(relatedness.link)
    group.net <- matrix(0,N,N)
    groups <- output.relat.samples[[which.net]][[which.sample]]$groups
    for (i in 1:length(groups)) {
      group.net[groups[[i]],groups[[i]]] <- sample(c(0,1),length(groups[[i]])^2,prob=c(0.9,0.1),replace=TRUE)
    }
    composite.network <- 10000*group.net + 10000*relatedness.link
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
    if(which.sample == 1){ 
      title(ylab="(a) Pedigree", font=2)
      #axis(2,labels=F, line=1)
    }  
    
    
    
    ## Plot relatedness network
    graph3 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.network,mode="undirected",diag=FALSE,weighted=TRUE)
    
    # Define coordinates for plotting group members nodes closer
    relatedness.network <- output.relat.samples[[which.net]][[which.sample]]$relatedness.network
    N <- nrow(relatedness.network)
    group.net <- matrix(0,N,N)
    groups <- output.relat.samples[[which.net]][[which.sample]]$groups
    for (i in 1:length(groups)) {
      group.net[groups[[i]],groups[[i]]] <- sample(c(0,1),length(groups[[i]])^2,prob=c(0.9,0.1),replace=TRUE)
    }
    composite.network <- 10000*group.net + 10000*relatedness.network
    graph.composite.network <- graph.adjacency(composite.network,mode="undirected",diag=FALSE,weighted=TRUE)
    coords <- layout_with_fr(graph.composite.network,weights=E(graph.composite.network)$weight)
    
    plot(graph3,
         edge.width=10*(exp(E(graph3)$weight)-1)^2,
         vertex.size=(as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1)*3,
         vertex.label="",         
         vertex.shape="fcircle", 
         vertex.frame.color="black",
         vertex.frame.width=0.5,
         vertex.color=c("black","red")[as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1],
         mark.groups=output.relat.samples[[which.net]][[which.sample]]$groups,
         layout=coords) 
    if(which.sample == 1){ 
      title(ylab="(b) Relatedness", font=2)
      #axis(2,labels=F, line=1)
    }
    
    
    
    ## Plot social network, using relatedness network to define clusters
    graph8 <- graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.network,mode="undirected",diag=FALSE,weighted=TRUE)
#    modul <- cluster_walktrap(graph8, weights = E(graph8)$weight, steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE)
    modul <- cluster_fast_greedy(graph8, weights = E(graph8)$weight, merges = TRUE, modularity = TRUE, membership = TRUE)
    
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
    composite.network <- 1000000*group.net + 1000*social.network
    graph.composite.network <- graph.adjacency(composite.network,mode="undirected",diag=FALSE,weighted=TRUE)
    coords <- layout_nicely(graph.composite.network,weights=E(graph.composite.network)$weight)
    
    # filtering network edges for better visualization
    edge.filter = params[which.net,4]
    coop.total.tmp2 <- coop.total.tmp
    coop.total.tmp2[coop.total.tmp2<edge.filter] <- 0
    
    # plot social network at the given time step with clustering from relatedness (colors)
    graph9 <- graph.adjacency(coop.total.tmp2,mode="undirected",diag=FALSE,weighted=TRUE)
    
    #V(graph9)$color <- modul$membership + 1

    #com=modul$membership
    ##comcolor <- rainbow(n=length(groups(modul)), alpha=1, start=1/6, end=4/6)
    ##comcolor <- rainbow(length(levels(as.factor(com))))[com]
    ##comcolor <- blue2green(length(groups(modul)))
    ##comcolor <- colorRampPalette(c(rgb(1,0,0,1), rgb(1,0,0,0)), alpha = TRUE)(length(groups(modul)))
    #comcolor <- colorRampPalette(c("white", "red"), bias=10)( length(groups(modul)) )
    #V(graph9)$color <- comcolor[com]
    V(graph9)$color <- colorRampPalette(c("white", "yellow", "red"), bias=10)( length(levels(as.factor(modul$membership))) )[modul$membership]
    
    plot(graph9,
   #     mark.groups=output.relat.samples[[which.net]][[which.sample]]$groups,
         edge.width=((E(graph9)$weight)),
         vertex.shape="fcircle", 
         vertex.frame.color="black",
         vertex.frame.width=0.5,
         vertex.size=(as.numeric(rowSums(output.relat.samples[[which.net]][[which.sample]]$coop.NOW)>=1)+1)*5,
         vertex.label="",
         edge.color="grey",
         layout=coords)
    #axis(1,labels=F, line=1)
    
    if(which.sample == 1){ 
      title(ylab="(c) Social", font=2)
      #axis(2,labels=F, line=1)
    }
    #mtext(text="Time step",side=1,line=1,outer=TRUE)
  } 

