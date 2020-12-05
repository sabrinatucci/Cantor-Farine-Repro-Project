# Cantor-Farine-Repro-Project

"Simple foraging rules in competitive environments can generate socially structured populations" (Cantor & Farine 2018)

Project goal: Reproduce Figure 3: Evolution of pedigree, relatedness, and social relationships in the reproductive model. 

12/4 UPDATE
-Successfully produced figure 3, working on progress report 

11/22 UPDATE 
- Combing through all of the different script files to determine which ones I will need and which ones I won't
- Re-did my github repository to fix reoccurring problems 
- Successfully ran the code scripts "set-up", "models" and "functions" (it took a long time for it to run)
- I ran the code for Figure 3, but got "Error in graph.adjacency(output.relat.samples[[which.net]][[which.sample]]$relatedness.link,  : 
  object 'output.relat.samples' not found". So I went back a ran the piece of code that had that line and was able to get it to produce a figure BUT I can't see the figure. It won't let me see it in R, and when I try to open it on my desktop it says that the file cannot be open because it may be damaged or use a file format that preview doesn't recognize 
  
  TO DO NEXT
  - Figure out how to get the code to successfully produce the figure 









General Info
  - Simulations were run for the parameter space area representing small population and resource patch sizes (see Figures 2b and 4b). The pedigree networks (a), showing all individuals that are alive at a given time step (x-axis), show that related individuals (nodes connected by links) are frequently found in foraging groups. Red nodes (and shading) indicate individuals currently part of the single foraging group. Similarly, the relatedness networks (b) show that individuals within social groups are often highly related (the thicknesses of links are proportional to their relatedness). However, the network depicting individuals are connected by the proportion of times they have been part of a specialized foraging group (c) suggests that foraging groups often contain individuals from different genetic lineages (here node color represents unique genetic lineages; for better visualization links whose weights <0.3 are filtered out). The simulations are based on a population size N = 40, resource patch size R = 15, initial connectivity T = 0.2, and run for 1,000 time steps, and these patterns are consistent across the parameter space (Figure S3)
  
  Paper Overview
  - Paper goal: simulate the emergence of foraging groups from simple interaction rules using an agent-based modeling framework 
  - Investigate whether the propensity to forage with the same conspecifics if the last attempt was successful is sufficient to produce (a) populations exhibiting resource-use specialization, and (b) structured patterns of relatedness among group members. 
  - Built two agent-based models, a baseline model and a reproductive model, in R 3.2.0 (R Development Core Team, 2017) 
  - Packages used: igraph, sna (Butts, asnipe and locfit (Loader, 2007) packages. 
  - All models, codes, and data are available at the online repository https://bitbucket.org/maucantor/coopgroup/src.
  
  Reasons for picking this paper:
  - Uses the same packages (igraph, asnipe) as the last figure I reproduced but adds complexity with more variables, more descriptive relationships and an aesthetic element.
  - Depicts social networks in unique ways 
  - This paper has 3 other interesting figures that I would like to reproduce if time permits. 

Next Steps
- Clone the bitbucket repository into github --> is this possible? or will I have to independetly download them?
- See if the code works 
