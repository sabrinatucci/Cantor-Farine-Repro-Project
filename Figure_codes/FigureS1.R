# Figure S1: 3d plots for model1 across parameter space

# Load packages and data
load(paste(getwd(), "/data/2_sensitivity_simulation.RData", sep=""))
source("setup.R")


# Plotting the 4 corners of the parameter space represented in the surface plot, Fig2

par(mfrow=c(4,3), mar=c(1,1,1,1))

# Area (b) in Fig2 b: small population, small resource
# N=40, r=5, tprob=0.2
z40.5n = z.plot(input.data=output.sensi, variable="n.groups", plot.parameters=c(40,5,0.2), x.limit=10, y.limits=c(1,100))
z40.5g = z.plot(input.data=output.sensi, variable="group.size", plot.parameters=c(40,5,0.2), x.limit=20, y.limits=c(1,100))
z40.5p = z.plot(input.data=output.sensi, variable="payoffs", plot.parameters=c(40,5,0.2), x.limit=NA, y.limits=c(1,100))

# Area (c) in Fig2 c: small population, large resource
# N=40, r=35, tprob=0.2
z40.35n = z.plot(input.data=output.sensi, variable="n.groups", plot.parameters=c(40,35,0.2), x.limit=10, y.limits=c(1,100))
z40.35g = z.plot(input.data=output.sensi, variable="group.size", plot.parameters=c(40,35,0.2), x.limit=20, y.limits=c(1,100))
z40.35p = z.plot(input.data=output.sensi, variable="payoffs", plot.parameters=c(40,35,0.2), x.limit=10, y.limits=c(1,100))

# Area (d) in Fig2 d: large population, large resource
# N=120, r=35, tprob=0.2
z120.35n = z.plot(input.data=output.sensi, variable="n.groups", plot.parameters=c(120,35,0.2), x.limit=10, y.limits=c(1,100))
z120.35g = z.plot(input.data=output.sensi, variable="group.size", plot.parameters=c(120,35,0.2), x.limit=20, y.limits=c(1,100))
z120.35p = z.plot(input.data=output.sensi, variable="payoffs", plot.parameters=c(120,35,0.2), x.limit=NA, y.limits=c(1,100))

# Area (e) in Fig2 e: large population, small resource
# N=120, r=5, tprob=0.2
z120.5n = z.plot(input.data=output.sensi, variable="n.groups", plot.parameters=c(120,5,0.2), x.limit=10, y.limits=c(1,100))
z120.5g = z.plot(input.data=output.sensi, variable="group.size", plot.parameters=c(120,5,0.2), x.limit=20, y.limits=c(1,100))
z120.5p = z.plot(input.data=output.sensi, variable="payoffs", plot.parameters=c(120,5,0.2), x.limit=10, y.limits=c(1,100))




# Alternative plot: other parameter combinations

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

