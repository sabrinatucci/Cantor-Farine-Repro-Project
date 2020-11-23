#' Model 1: foraging group evolution
#' 
#' @param N Integer. Population size: number of nodes in the social network
#' @param resource.size Integer. Size of resource patch
#' @param n.reps Integer. Time steps for the model run (number of repetition)
#' @param tprob Initial density (connectance) of the network. Proportion of realized links in relation to all possible links. Links will be set randomly among individuals
#' @param type string: Type of resource allocation to individuals in groups. If \code{type=size-based} it calculates individual payoffs in which the group-level resource share is determined by the size of the groups, given by \code{resource.size*(group.lengths^2)/sum(group.lengths^2)},  considering that larger groups outcompete smaller ones; if \code{type='equal'}, it considers no disproportionate individual payoffs, that is, there is equal allocation of resource to all individuals irrespective of their group size, as in \code{resource.size*group.lengths/sum(group.lengths)}. See details in the help file of the \code{calculate.payoffs()} function.
#' @param output "model", "exclusivity", "sensitivity"
#' @return If \code{output == "model"}, the function returns a list with 5 objects: 
#'    \code{coop.total}: matrix with all interactions across time steps to be plotted as a network; 
#'    \code{mean.payoffs}: vector with average individual payoff per time step; 
#'    \code{n.groups}: vector with number of groups across per time steps 
#'    \code{mean.roup.size}: vector with the  average group size per time step; 
#'    \code{exclusivity}: vector with exclusivity of the groups per time step;
#'    \code{groups}: list with all foraging groups at the last time step; 
#'    \code{coop.NOW}: matrix with social interactions among all alive individuals at the last time step. 
#' 
#' If \code{output == "exclusivity"}, the function returns a matrix with 4 columns ("N","r","tprob","Exclusivity") to be used in the surface plot. 
#' 
#' If \code{output == "sensitivity"}, the function returns a matrix with 8 columns("time.step", "N","r","tprob","Exclusivity", "Mean.payoffs","n.groups", "Mean.group.size") to be used in the sensitivity analysis
#' @export
model1 <- function(N, resource.size, n.reps, tprob, type="size-based", output="model"){
  
  # Setting initial directed network: start with a random set of edges; no loops
  coop <- rgraph(N,tprob=tprob); diag(coop) <- NA
  # Network with average number of foraging association across runs
  coop.total <- matrix(0,nrow=N,ncol=N)  
  # Mean individual payoffs
  mean.payoff <- rep(NA,n.reps)
  # All formed groups
  groups.total <- list()
  # Exclusivity
  group.exclusivity <- rep(NA,n.reps)
  
  # Model run
  for (zz in 1:n.reps) {
    coop.NOW <- matrix(0,nrow=N,ncol=N)
    
    # extract pairs of cooperators (reciprocal links) out of the network
    ABind <- cooperators(coop,N)
    # form groups based on chain rule (A->B, B->C, then A->C)
    groups <- identify.groups(ABind)
    # calculate per capita payoff based on given resource patch size
    if (length(groups) >= 1) {
      payoffs <- calculate.payoffs(groups,resource.size,type) # calculate the group-level resource share based on group size
      # individual payoffs
      inds.payoffs <- rep(0,N)
      for (i in 1:length(groups)) {
        inds.payoffs[groups[[i]]] <- payoffs[i]
        a <- expand.grid(groups[[i]],groups[[i]])
        a <- a[which(a[,1] != a[,2]),]
        coop.NOW[as.matrix(a)] <- 1  
        coop.total[as.matrix(a)] <- coop.total[as.matrix(a)]+1
      }
    }
    
    #update network
    coop <- update.coop.network(coop,inds.payoffs)
    
    if(output != "exclusivity"){
      # Mean nonzero individual payoff
      mean.payoff[zz] <- mean(inds.payoffs[which(inds.payoffs!=0)])
      # Number and size of groups for each run
      groups.total[[zz]] <- groups
      # Exclusivity: proportion of times individuals were part of the foraging group for each run
      group.exclusivity[zz] <- exclusivity(initial.size=N, runs=zz, patch=resource.size, all.ties=coop.total)
    }
  }
  
  # OUTPUTS
  if(output != "exclusivity"){
    group.number <- groups.number(all.groups=groups.total, runs=n.reps)
    group.size <- average.group.size(groups.total)
    if(output == "model"){
      result.model <- list(coop.total, mean.payoff, group.number, group.size, group.exclusivity, coop.NOW, groups)
      names(result.model) <- c("coop.total", "mean.payoffs", "n.groups", "mean.group.size", "exclusivity", "coop.NOW", "groups")
      return(result.model)
    } else { # output == 'sensitivity'
      result.model <- cbind(1:n.reps, rep(N,n.reps), rep(resource.size,n.reps), rep(tprob,n.reps), group.exclusivity, mean.payoff, group.number, group.size)
      colnames(result.model) <- c("time.step", "N","r","tprob","Exclusivity", "Mean.payoffs","n.groups", "Mean.group.size")
      return(result.model)
    }
  } else { # output == "exclusivity"
    group.exclusivity <- exclusivity(initial.size=N, runs=(n.reps-N), patch=resource.size, all.ties=coop.total)
    result.model <- c(N, resource.size, tprob, group.exclusivity)
    return(result.model)
  }
}







#' Model 2: multi-generational foraging group evolution (adding relatedness)
#' 
#' @param N Integer. Population size: number of nodes in the social network
#' @param resource.size Integer. Size of resource patch
#' @param n.reps Integer. Time steps for the model run (number of repetition)
#' @param tprob Initial density (connectance) of the network. Proportion of realized links in relation to all possible links. Links will be set randomly among individuals
#' @param type string: Type of resource allocation to individuals in groups. If \code{type=size-based} it calculates individual payoffs in which the group-level resource share is determined by the size of the groups, given by \code{resource.size*(group.lengths^2)/sum(group.lengths^2)},  considering that larger groups outcompete smaller ones; if \code{type='equal'}, it considers no disproportionate individual payoffs, that is, there is equal allocation of resource to all individuals irrespective of their group size, as in \code{resource.size*group.lengths/sum(group.lengths)}. See details in the help file of the \code{calculate.payoffs()} function.
#' @param output If \code{output = "model"}, the model runs until the end and the results given regards the last time step; 
#'               If \code{output = "samples"}, the results are given for 5 time steps, equally spaced, during the simulation. E.g. if the total simulation time is 1000, there will be outputs for t = (200, 400, 600, 800, 1000).
#'               If \code{output = "relatedness"}, the model runs as in \code{output = "samples"} but also computes relatedness among foraging group members and non-members individuals
#'               If \code{output = "sensitivity"}, the model runs as in \code{output = "model"} but outputs the ratio between the mean relatedness among group members and the mean relatedness among non-members at the last time step
#'               
#' @return If \code{output == "model"}, the function returns the output of a complete model run in a list with 5 objects: 
#'    \code{coop.NOW}: matrix with social interactions among all alive individuals at the last time step; 
#'    \code{groups}: list with all foraging groups at the end of simulation; 
#'    \code{ids}: data frame with 4 columns: Individual identification (ID), its parent indentification (Parent), if it is alive (1) or not (0), its age in number of time steps; 
#'    \code{alive}: matrix with all alive individuals, to be used for plotting social network; 
#'    \code{coop.total}: matrix with all interactions across time steps to be plotted as a network; 
#'    
#'        If \code{output == "samples"}, the function returns a list with 5 objects, each of which is a list with the same outputs when \code{output == "model"} for 5 snapshots along the simulation
#'        
#'        If \code{output = "relatedness"}, the function returns the same outputs as when \code{output == 'samples'}, but with five additional objects:
#'    \code{relatedness.link}: matrix with binary relatedness links among all alive individuals at the last time step; 
#'    \code{relatedness.network}: matrix with relatedness interactions (based on geodesic distances) among all alive individuals at the last time step, which can be plotted as a social network; 
#'    \code{relatedness.link.total} matrix with sum of all binary relatedness links that existed among individuals during simulation; 
#'    \code{relatedness.network.total}: matrix with sum of all relatedness values during simulation,  which can be plotted as a social network; 
#'    \code{all.relatedness}: matrix with six columns: $mean.coop.relatedness: mean relatedness among individuals of the emergent cooperative groups per time step; $mean.noncoop.relatedness: mean relatedness among all non-cooperative individuals per time step; $mean.random.relatedness: mean random relatedness form permutations, using a randomisation test with 100 iterations; $low2.5CIrandom: 2.5% confidence interval of the random distribution of relatedness; $upper97.5CIrandom: 97.5% confidence interval of the random distribution of relatedness; $relat.ratio: ratio between the average relatedness among individuals of the foraging group, divided by the mean relatedness from the random simulated data
#'    
#'        If \code{output = "sensitivity"}, the function returns a matrix with 5 columns ("N","r","tprob","Relatedness.ratio", "logRelatedness.ratio") to be used in the sensitivity analysis. Relatedness ratio is calculated at the LAST time step only, and as the average relatedness among individuals of the emergent foraging groups divided by the average relatedness among all non-member individuals
#'        
#' @export
model2 <- function(N, resource.size, n.reps, tprob, type="size-based", output="model"){
  
  if(output == "samples" | output == "relatedness"){
    # defining 5 time steps across the model run to pull out results from
    SAMPLES <- seq(from=(n.reps/5), to=n.reps, by=n.reps/5)
    SAMPLES <- SAMPLES[order(SAMPLES)]
  } else {
    # output for the last time step only
    SAMPLES <- n.reps
  }
  
  # Setting initial directed network: start with a random set of edges; no loops
  coop <- rgraph(N,tprob=tprob)
  diag(coop) <- NA
  
  # Keep track of agents' information
  ids <- data.frame(ID=1:N,Parent=NA,Alive=1,Age=0)
  coop.total <- matrix(0,nrow=N,ncol=N)
  alive <- matrix(0,nrow=n.reps,ncol=N)
  samples.out <- list()
  permutation.relatedness <- matrix(NA,n.reps, 100) # rows=time steps, columns= 100 iterations for relatedness permutation test
  mean.permuted.relat <- list()
  all.relatedness <- data.frame(mean.coop.relatedness=numeric(), mean.noncoop.relatedness=numeric(),mean.random.relatedness=numeric(), low2.5CIrandom=numeric(), upper97.5CIrandom=numeric(), relat.ratio=numeric())
    
  
  # RUN MODEL 2
  for (zz in 1:n.reps) {
        
    # extract pairs of cooperators (reciprocal links) out of the network
    ABind <- cooperators(coop,N)
    # form groups based on chain rule (A->B, B->C, then A->C)
    groups <- identify.groups(ABind)
    # save cooperative network
    coop.NOW <- matrix(0,nrow=N,ncol=N)
    
    # calculate per capita payoff based on given resource patch size
    inds.payoffs <- rep(0,N)    
    if (length(groups) >= 1) {
      # individual payoff
      payoffs <- calculate.payoffs(groups,resource.size,type) # calculate the group-level resource share based on group size
      for (i in 1:length(groups)) {
        inds.payoffs[groups[[i]]] <- payoffs[i]
        a <- expand.grid(groups[[i]],groups[[i]])
        a <- a[which(a[,1] != a[,2]),]
        coop.NOW[as.matrix(a)] <- 1
        a[,1] <- ids$ID[ids$Alive==1][a[,1]]
        a[,2] <- ids$ID[ids$Alive==1][a[,2]]
        coop.total[as.matrix(a)] <- coop.total[as.matrix(a)]+1     
      }
    }
    # update cooperative network
    coop <- update.coop.network(coop,inds.payoffs)
        
    # update agents's age and who's alive   
    ids$Age[which(ids$Alive==1)] <- ids$Age[which(ids$Alive==1)] + 1
    alive[zz,which(ids$Alive==1)] <- 1
    
    # calculate relatedness among members of foraging groups and non-members
    if(output == "relatedness"){
        # relatedness among all agents
        rel.link.total <- get.relatedness.link(ids)
        rel.network.total <- get.relatedness.network(rel.link.total)
        diag(rel.network.total) <- NA # make sure to remove self relatedness for non-members     
        
        # relatedness among all alive agents
        rel.link.alive <- rel.link.total[which(ids$Alive==1),which(ids$Alive==1)] 
        rel.network.alive <- get.relatedness.network(rel.link.alive)
        diag(rel.network.alive) <- NA    
              
        # average relatedness among foraging group members, among non-members, and at random
          # if there are no groups, skip this part
          if(length(groups)==0){ 
            all.relatedness[zz,] <- NA
          } else {
          # get the individuals for each group, and the relatedness among them
            aux.rel <- list()        
            for (i in 1:length(groups)) {
              coops=unlist(groups[i])
              aux=rel.network.total[coops,coops]
              aux.rel[[i]]=aux[which(lower.tri(aux))]
            }
            # average relatedness among foraging group members
            all.relatedness[zz,1] <- mean(unlist(aux.rel),na.rm=TRUE)
            # average relatedness among all non-members
            allcoops=unlist(groups)          
            all.relatedness[zz,2] <- mean(rel.network.alive[-allcoops, -allcoops], na.rm=T) 
    
            # Permutation test: calculating relatedness expected by chance to compare with foraging group members
            for(k in 1:ncol(permutation.relatedness)){             
              # form groups of the same size of a given foraging group by resampling alive individuals from the population
                group.sample = resample.groups(groups = groups, inds = ids[which(ids$Alive==1), ]$ID)     
                for(g in 1:length(group.sample)){
                  mean.permuted.relat[[g]] = mean(rel.network.total[ group.sample[[g]], group.sample[[g]] ], na.rm=T)            
                }
                # extract the observed total relatedness among the individuals of the random group and save their mean relatedness
                permutation.relatedness[zz,k] = mean(unlist(mean.permuted.relat), na.rm=T)
                mean.permuted.relat <- list()
            }               
            # save mean random relatedness and the 95% CI 
            all.relatedness[zz,3] <- mean(permutation.relatedness[zz,])
            all.relatedness[zz,4:5] <- as.numeric(quantile(permutation.relatedness[zz,], probs=c(0.025, 0.975), na.rm=TRUE, type=2))
          
            # calculate the ratio of observed relatedness of the real foraging group to the mean random relatedness
            # (mean relatedness among all the real foraging groups at a given time step) DIVIDED BY (mean random relatedness at a given time step)
            all.relatedness[zz,6] <- ((all.relatedness[zz,1]) / (all.relatedness[zz,3]))        
          } 
        }
    

    # Reproduction (stochastic)
    repro <- sapply(reproduce.die(ids$Age[which(ids$Alive==1)]),function(x) { sample(c(0,1),1,prob=c(1-x,x))})
    
      # Update agents post reproduction
      if (sum(repro) >= 1) {
        new.ids <- data.frame(ID=(max(ids$ID)+1):((max(ids$ID)+sum(repro))),Parent=ids$ID[which(ids$Alive==1)][which(repro==1)],Alive=1,Age=0)
        ids <- rbind(ids,new.ids)
        new.rows <- coop[which(repro==1),]
        new.rows[is.na(new.rows)] <- 0
        coop <- rbind(coop,new.rows)
        new.cols <- coop[,which(repro==1)]
        new.cols[is.na(new.cols)] <- 1
        coop <- cbind(coop,new.cols)
        diag(coop) <- NA
        alive <- cbind(alive,matrix(0,nrow=n.reps,ncol=sum(repro)))
        coop.total <- rbind(coop.total,matrix(0,nrow=sum(repro),ncol=ncol(coop.total)))
        coop.total <- cbind(coop.total,matrix(0,ncol=sum(repro),nrow=nrow(coop.total)))
      }
    
    # Death (stochastic)
    die <- sapply(reproduce.die(ids$Age[which(ids$Alive==1)]),function(x) { sample(c(0,1),1,prob=c(1-x,x))})
      
      # Update agents post death
      if (sum(die) >= 1) {
        coop <- coop[which(die!=1),which(die!=1)]
        ids$Alive[which(ids$Alive==1)[which(die==1)]] <- 0
      }
    
    # Update population size
    N <- sum(ids$Alive)
    
    # Save samples output   
    if(output=="model" | output=="samples"){
      if (zz %in% SAMPLES) {
        samples.out[[which(SAMPLES==zz)]] <- list(
          coop.NOW=coop.NOW,
          groups=groups,
          ids=ids,
          alive=alive,
          coop.total=coop.total)
      }
    }
    
    if(output=="relatedness"){
      if (zz %in% SAMPLES) {
        samples.out[[which(SAMPLES==zz)]] <- list(
          coop.NOW=coop.NOW,
          groups=groups,
          ids=ids,
          alive=alive,
          coop.total=coop.total,
          relatedness.link=rel.link.alive, 
          relatedness.network=rel.network.alive, 
          relatedness.link.total=rel.link.total, 
          relatedness.network.total=rel.network.total, 
          all.relatedness=all.relatedness)
      }
    }
  }
  
  # OUTPUT
  if(output == "sensitivity"){
    # Calculate relatedness ratio at the end of simulation for alive individuals only
    rel.link.total <- get.relatedness.link(ids)
    rel.network.total <- get.relatedness.network(rel.link.total)
    rel.link.alive <- rel.link.total[which(ids$Alive==1),which(ids$Alive==1)] 
    rel.network.alive <- get.relatedness.network(rel.link.alive)
    diag(rel.network.alive) <- NA    
    
    # get individuals for each group, calculate their relatedness, then average relatedness among all groups
    if(length(groups)==0){
      # if there are no groups, skip 
      relat.ratio <- NA
    } else {
      # get the individuals for each group, and the relatedness among them
      aux.rel <- list()
      for (i in 1:length(groups)) {
        coops=unlist(groups[i])
        aux=rel.network.total[coops,coops]
        aux.rel[[i]] <- aux[which(lower.tri(aux))]
      }
      allcoops = as.numeric(unlist(groups))
      # calculate the ratio relatedness among all foraging group members individuals to the relatedness among all alive non-members
      relat.ratio <- (mean(unlist(aux.rel),na.rm=TRUE)) /  (mean(rel.network.alive[-allcoops, -allcoops], na.rm=TRUE))
    } 
    # outputs data for surface plot
    result.model <- c(N, resource.size, tprob, relat.ratio, log(relat.ratio))
    return(result.model)
    
  } else {
    # output the list with 5 or 10 object
    for(i in 1:length(SAMPLES)){ names(samples.out)[i] <- paste("Time.step.", SAMPLES[i], sep="")}
    return(samples.out)  
  }
}

