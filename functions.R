#' Returns pairs of cooperators
#' @param coop matrix: square network adjacency matrix in which aij =1 when individual in row i interacts with individual in column j, and aij=0 otherwise
#' @param N scalar: number of individuals in the network
#' @return ABind data.frame: pairs of individuals with reciprocal links in the \code{coop} network
cooperators <- function(coop, N) {
	ABind <- expand.grid(1:N,1:N)
	ABind <- ABind[which(ABind[,1] < ABind[,2]),]
	ABind <- ABind[order(ABind[,1]),]
	AB <- t(coop)[lower.tri(coop)]
	BA <- coop[lower.tri(coop)]
	return(ABind[which(AB == 1 & BA == 1),])
}





#' Forms groups based on chain rule (A->B, B->C, then A->C)
#'@param ABind data.frame: pairs of individuals with reciprocal links in the \code{coop} network
#'@return groups list: vectors with individuals who are members of foraging groups
identify.groups <- function(ABind) {
	groups <- list()
	g <- 1
	inds <- unique(c(ABind[,1],ABind[,2]))
	
	while (length(inds) >= 1) {
		grp.tmp <- inds[1]
		a <- which(ABind[,1] %in% grp.tmp | ABind[,2] %in% grp.tmp)

		while (length(a) >= 1) {
			grp.tmp <- c(grp.tmp, unique(c(ABind[a,1],ABind[a,2])))
			grp.tmp <- unique(grp.tmp)
			inds <- inds[which(!(inds %in% grp.tmp))]
			ABind <- ABind[-a,]
			a <- which(ABind[,1] %in% grp.tmp | ABind[,2] %in% grp.tmp)
		}
		
		groups[[g]] <- grp.tmp
		g <- g + 1
	}
	return (groups)
}





#' Calculates payoff for individuals in all groups based on given resource patch size
#' @param groups list: list of vectors with individuals who are members of foraging groups
#' @param resource.size scalar: size of the resource patch
#' @param type string: Type of resource allocation to individuals in groups. If \code{type='size-based'} it calculates individual payoffs in which the group-level resource share is determined by the size of the groups, given by \code{resource.size*(group.lengths^2)/sum(group.lengths^2)},  considering that larger groups outcompete smaller ones; if \code{type='equal'}, it considers no disproportionate individual payoffs, that is, there is equal allocation of resource to all individuals irrespective of their group size, as in \code{resource.size*group.lengths/sum(group.lengths)}.
#' @retunr per.capita.payoff number: payoff for individuals in the group 
calculate.payoffs <- function(groups, resource.size, type="size-based") {
	group.lengths <- unlist(lapply(groups,length))
	if (type == "size-based") {
		group.payoffs <- resource.size*(group.lengths^2)/sum(group.lengths^2)
	} else { # e.g. type=="equal"
		group.payoffs <- resource.size*group.lengths/sum(group.lengths)
	}
	per.capita.payoff <- group.payoffs / group.lengths
	return(per.capita.payoff)
}	





#' Updates network ties based on individual payoffs on foraging in group on a resource patch 
#' @param coop matrix: square network adjacency matrix in which aij =1 when individual in row i interacts with individual in column j, and aij=0 otherwise
#' @param inds.payoffs vector: payoffs for all individuals in the groups of current time step t, as calculated by \code{per.capita.payoff}
#' @param p.stochastic number: probability of a random tie to be added to the network
#' @return coop matrix: updated network adjacency matrix with added and removed ties according to the individuals payoffs
update.coop.network <- function(coop,inds.payoffs,p.stochastic=0.0001) {
	N <- length(inds.payoffs)
	for (i in which(inds.payoffs > 0)) {
		# reduce cooperators
		if (inds.payoffs[i] < 1) {
			if (sum(coop[i,],na.rm=TRUE) > 0) {
				if (sum(coop[i,],na.rm=TRUE) > 1) {
					coop[i,sample(which(coop[i,]==1),1)] <- sample(c(0,1),1,prob=c(0.8,0.2))
				} else {
					coop[i,which(coop[i,]==1)] <-  sample(c(0,1),1,prob=c(0.5,0.5))
				}
			}
		}
		if (inds.payoffs[i] > 1) {
			if (sum(coop[i,],na.rm=TRUE) < (N-1)) {
				if (sum(coop[i,]==0,na.rm=TRUE) > 1) {
					coop[i,sample(which(coop[i,]==0),1)] <- sample(c(1,0),1,prob=c(0.1,0.9))
				} else {
					coop[i,which(coop[i,]==0)] <-  sample(c(1,0),1,prob=c(0.1,0.9))
				}
			}
		}
	}

	# Stochastic addition of edges
	ABs <- expand.grid(1:N,1:N)
	ABs <- ABs[which(ABs[,1] != ABs[,2]),]
	new <- ABs[sample(nrow(ABs),round((N^2)*p.stochastic)),]
	coop[new[,1],new[,2]] <- 1
	diag(coop) <- NA

	# Stochastic loss of edges
	#ABs <- expand.grid(1:N,1:N)
	#ABs <- ABs[which(ABs[,1] != ABs[,2]),]
	#new <- ABs[sample(nrow(ABs),round((N^2)*p.stochastic)),]
	#coop[new[,1],new[,2]] <- 0

	return (coop)
}





#' Calculates Exclusivity, as the proportion of times individuals were part of the foraging group
#' @param initial.size scalar: network size = \code{N}
#' @param runs scalar: number of model runs = \code{n.reps}
#' @param patch scalar: size of the resource patch =  \code{resource.size}
#' @param all.ties matrix: Network with average number of times individuals were part of the same emergent foraging group across runs = \code{coop.total}
#' @return exclusivity number
exclusivity <- function(initial.size=N, runs=n.reps, patch=resource.size, all.ties=coop.total){
  coop.total.tmp <- all.ties/runs
  degs <- rowSums(coop.total.tmp)
  ids <- 1:initial.size
  ids <- ids[rev(order(degs))]
  sum.top.r <- sum(coop.total.tmp[ids[1:patch],ids[1:patch]],na.rm=T)
  exclusivity <- sum.top.r/sum(coop.total.tmp,na.rm=T)
  return(exclusivity)
}





#' Calculates number of groups across runs
#' @param groups.total list: each element contain all groups (as lists) in each model run
#' @param runs scalar: number of model runs = \code{n.reps}
#' @return summary vector: number of emergent foraging groups in each model run
groups.number <- function(all.groups=groups.total, runs=n.reps){  
  summary <- numeric()
  for(zz in 1:runs){
    summary[zz] <- length(all.groups[[zz]])
  }
  return(summary)
}





#' Calculates average group size per model run
#' @param groups.total list: each element contain all groups (as lists) in each model run
#' @return avg vector: average group size per model run
average.group.size <- function(groups.total) {
  gs <- length(groups.total)
  avg <- rep(NA,gs)
  
  for (i in 1:gs) {
    avg[i] <- sum(unlist(lapply(groups.total[[i]],length)))/length(groups.total[[i]])
    
  }
  return(avg)
}





#' Calculates probability of group sizes during simulation
#' @param groups list: list of vectors with individuals who are members of foraging groups
#' @param N scalar: number of individuals in the network
#' @return group.probs vector: probability of group sizes during simulation for building 3D plots, see \code{z.plot}
get.group.probability <- function(groups, N) {
  group.probs <- rep(0,N+1)
  groups.table <- table(groups)
  group.probs[as.numeric(names(groups.table))+1] <- groups.table
  group.probs <- group.probs/sum(is.finite(groups))
  return(group.probs)
}





#' Calculates probability of individual payoffs during simulation
#' @param groups list: list of vectors with individuals who are members of foraging groups
#' @param bins vector: size of histogram bins from which the probability will be calculated
#' @return group.probs vector: probability of individual payoffs during simulation for building 3D plots, see \code{z.plot}
get.payoff.probability <- function(groups, bins) {
  a <- hist(groups,breaks=bins,plot=FALSE,include.lowest=TRUE)
  group.probs <- a$counts/sum(a$counts)
  return(group.probs)
}





#' Calculates probability of an agent to reproduce or die
#' @param age vector: age of individuals in time steps
#' @return Probability of an agent to reproduce or die
reproduce.die <- function(age) {
  return(0.001/(1+exp(-0.1*age+5)))
}





#' Records the relatedness between agents
#' @param ids data.frame: each row contain an individual, and the columns gives its ID, Parent ID, if it's Alive or not at the current time t, and its Age in time steps
#' @return relatedness.link matrix: square and binary matrix giving pedigree of individuals based on the parent-offspring relationships
get.relatedness.link<- function(ids) {
  N <- nrow(ids)
  relatedness.link <- matrix(0,nrow=N,ncol=N)
  
  for (i in 1:N) {
    offspring <- which(ids$Parent == i)
    if (length(offspring) > 0) {
      relatedness.link[i,offspring] <- 1
      relatedness.link[offspring,i] <- 1
    }
  }
  return(relatedness.link)
}





#' Creates adjacency matrix for relatedness network based on geodesic paths
#' @param relatedness.link matrix: square and binary matrix giving pedigree of individuals based on the parent-offspring relationships
#' @return relatedness.network matrix: square and weighted matrix giving relatedness among all individuals based on geodesic paths. Tie weights are calculated based on the pedigree relationships using lij = 1/2^(dij), where lij is the relatedness between individuals i and j, and dij is the shorted path length connecting the agents i and j (but lij is set to 0 if they are not connected). 
get.relatedness.network <- function(relatedness.link) {
  relatedness.network <- 1/(2^geodist(relatedness.link)$gdist)
  return(relatedness.network)
}





#' Running a given model across a given parameter space
#' @param inputs vector: give the four main input parameters for the models: N = population size; r = resource patch size, tprob = Initial density (connectance) of the network; type = type of group resource share. Note the parameter 'n.reps' will be automatically: if you run the \code{model="model1"} then \code{n.reps = N*5}; if \code{model="model2"} then \code{n.reps = 1000}. See help files of \code{model1, model 2} for details on the parameters.
#' @param model string: "model1" to run \code{model1()}; "model2" to run \code{model2}
#' @param output string: for \code{model1()} the possible outputs are c("model", "exclusivity", "sensitivity"); for \code{model2()}, the output options are c("model", "samples", "relatedness", "sensitivity"); see \code{model1, model 2} 
#' @details see \code{model1, model 2} 
#' @return returns the model output as defined in the \code{output} paramenter
simulation <- function(inputs, model, output){
  cat("running parameters:", as.numeric(inputs),"\n")

  # assign input parameters
  N <- as.numeric(inputs[1])
  resource.size <- as.numeric(inputs[2])
  tprob <- as.numeric(inputs[3])
  type <- as.character(inputs[,4])
  n.reps <- 5*N
  
  
  #  # when resource patch size is <= population size
  #  if (resource.size <= N) { 
  if(model=="model1") { 
    result <- model1(N, resource.size, n.reps, tprob, type, output)
  } else {
    n.reps <- 1000
    result <- model2(N, resource.size, n.reps, tprob, type, output)}
  #}
  return(result)
}





#' Internal function for 3D plot
trans3d <- function(x,y,z, pmat) {
  tr <- cbind(x,y,z,1) %*% pmat
  list(x = tr[,1]/tr[,4], y= tr[,2]/tr[,4])
}





#' Internal function for transparent 3D plot
add.alpha <- function(col, alpha=1){
if(missing(col))
  stop("Please provide a vector of colours.")
apply(sapply(col, col2rgb)/255, 2, 
      function(x) 
        rgb(x[1], x[2], x[3], alpha=alpha))  
}





#' Plots 3D graphs
#' @param input.data Output matrix from \code{model1, output == "sensitivity"}, which has 8 columns ("time.step", "N","r","tprob","Exclusivity", "Mean.payoffs","n.groups", "Mean.group.size") 
#' @param plot.parameters Vector defining the parameters to be plot in this order: N, r, tprob. For example \code{plot.parameters = c(40, 5, 0.2)} extracts from {input.data} all values for a network size N=40, resource patch size r=5, and initial connectanc tprob=0.2
#' @param variable. Which variable should be plotted? If \code{variable = n.groups} the z-axis
#' @param x.limit Integer. Upper limit of x-axis to be plot. For example, when plotting number of groups, \code{x.limit = 10} plot maximum of 10 groups
#' @param y.limits Vector defining the interval of time from the simulation (1000 time steps) to be plotted. For example \code{y.limits = c(1,50)} will plot the simulated data form time step 1 until time step 50.
#' @return z matrix with data to plot
#' @return 3D plot
z.plot <- function(input.data, plot.parameters, variable, x.limit, y.limits){
  
  plotcol="black"
  # define parameters
  N <- plot.parameters[1]
  r <- plot.parameters[2]
  tprob <- plot.parameters[3]
  N.max <- x.limit
  t.min <- y.limits[1]
  t.max <- y.limits[2]
  
  # prepare data
  input.data <- as.data.frame(input.data)
  input.data$ID <- 1:nrow(input.data)
  
  # define bins
  t.bins <- seq(t.min,t.max,1)
  if(variable != "payoffs") { N.groups.bins <- 0:N 
  } else { 
    N.groups.bins <- seq(from=0, to=max(input.data$Mean.payoffs)+0.5, 0.5) 
  }
  
  # defines z-axis based on given arguments
  z <- matrix(0,nrow=length(N.groups.bins),ncol=length(t.bins),byrow=FALSE)
  for (i in t.bins) {
    if(variable == "n.groups") { 
      label="Number of groups"
      input <- input.data[which(input.data$tprob==tprob & input.data$N==N & input.data$r == r & input.data$time.step==i),]$n.groups
      z[,i] <- get.group.probability(input,N)
      N.groups.bins <- 1:N.max
    }
    if(variable == "group.size") { 
      label="Average group size"
      input <- input.data[which(input.data$tprob==tprob & input.data$N==N & input.data$r == r & input.data$time.step==i),]$Mean.group.size
      z[,i] <- get.group.probability(round(input),N)
      N.groups.bins <- 1:N.max
    }
    if(variable == "payoffs") { 
      label="Average individual payoff"
      input <- input.data[which(input.data$tprob==tprob & input.data$N==N & input.data$r == r & input.data$time.step==i),]$Mean.payoffs
      z[,i] <- get.payoff.probability(input,c(N.groups.bins,max(N.groups.bins)+1)-0.25)
    }
  }
  
  # Define plot dimensions and colors
  if(variable != "payoffs") { 
    z2 <- z[1:N.max,] 
  } else { z2 <- z 
  }
  zmax <- max(z2)
  zmin <- 0
  minz <- min(z2,na.rm=T)
  nrz <- nrow(z2)
  ncz <- ncol(z2)
  nbcol <- 100
  jet.colors <- blue2green2red(nbcol)
  zfacet <- (z2[-1, -1] + z2[-1, -ncz] + z2[-nrz, -1] + z2[-nrz, -ncz])/4
  facetcol <- cut(zfacet,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
  zcol <- cut(z2,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
  
  # 3D plot
  res <- persp(N.groups.bins,t.bins, z2, theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors[facetcol], ylab="Time",xlab=label,zlab="Probability",zlim=c(zmin,zmax),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75)
  
  return(z)
}





#' Creates random groups by resampling alive individuals from the population
#' @param groups list: list of vectors with individuals who are members of foraging groups
#' @param inds vector: alive individuals from the population at time step t
#' @return groups.out list: list with vectors a random group of same size as \code{groups}
resample.groups <- function(groups, inds) {
  groups.out <- list()
  for (i in 1:length(groups)) {
    groups.out[[i]] <- sample(inds,length(groups[[i]]))
    inds <- inds[which(!(inds %in% groups.out[[i]]))]
  }
  return(groups.out)
}





#' Plots scatter plot with whiskers
#' @param x vector: values for x-axis
#' @param y vector: values for y-axis
#' @param yplus vector: upper whisker values
#' @param yminus vector: lower whisker values
#' @param ... arguments to be passed for \code{plot} function
#' @author Originally Charles Geyer, U.Chicago, early 1991; then Martin Mächler; package {sfsmisc}
errbar <- function (x, y, yplus, yminus, cap = 0.015, ylim = range(y, yplus, yminus), xlab = deparse(substitute(x)), ylab = deparse(substitute(y)), colwhisker='#0000001A',...) {
  plot(x, y, ylim = ylim, xlab = xlab, ylab = ylab, las=1,...)
  xcoord <- par()$usr[1:2]
  segments(x, yminus, x, yplus, col=colwhisker)
  smidge <- cap * (xcoord[2] - xcoord[1])/2
  segments(x - smidge, yminus, x + smidge, yminus, col=colwhisker)
  segments(x - smidge, yplus, x + smidge, yplus, col=colwhisker)
}





#' Plotting surface relatedness ratio
#' @param input data.frame: raw data for surface plot, usually with columns N,  r, tprob, Relatedness.ratio, logRelatedness.ratio,  ID
#' @param surf locfit object: surface function on \code{input} data using \code{locfit} function in {locfit} R package
#' @param xaxis string: label for x axis
#' @param yaxis string: label for y axis
#' @param zaxis string: label for z axis
plotsurf <- function(input, surf, yaxis, zaxis, xaxis){
  r.bins <- seq(min(input$r),max(input$r),1)
  N.bins <- seq(min(input$N),max(input$N),1)
  plotcol="black"
  z <- matrix(predict(surf,newdata=expand.grid(r=r.bins,N=N.bins),type="response"),nrow=length(r.bins),ncol=length(N.bins),byrow=FALSE)
  N.mat <- matrix(rep(N.bins,each=length(r.bins)),ncol=length(N.bins),nrow=length(r.bins))
  r.mat <- matrix(rep(r.bins,each=length(N.bins)),ncol=length(N.bins),nrow=length(r.bins),byrow=TRUE)
  zmax <- max(z)
  zmin <- min(z)
  minz <- min(z,na.rm=T)
  nrz <- nrow(z)
  ncz <- ncol(z)
  nbcol <- 100
  jet.colors <- blue2green2red(nbcol)
  jet.colors2 <- add.alpha(jet.colors,alpha=0.6)
  zfacet <- (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
  facetcol <- cut(zfacet,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
  zcol <- cut(z,breaks=seq(zmin,zmax,(zmax-zmin)/(nbcol)),labels=c(1:nbcol))
  res <- persp(r.bins, N.bins, z, theta=120, phi=30, shade=0.2, ticktype="detailed", expand=0.8, col=jet.colors2[facetcol], ylab=yaxis,xlab=xaxis,zlab=zaxis,zlim=c(zmin,zmax),nticks=5, border="black",lwd=0.1, cex.lab=1.2,ltheta = 235, lphi = 75)
}





#' Add parameter to igraph plot to change node border width
#' @author Gabor Csardi
#' @references http://lists.gnu.org/archive/html/igraph-help/2013-03/msg00030.html
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd, circles=size, add=TRUE, inches=FALSE)
         })
}
add.vertex.shape("fcircle", clip=igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1,vertex.frame.width=1))