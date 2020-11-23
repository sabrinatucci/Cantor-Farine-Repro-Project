# Initial setup

# Loading packages and functions
if(!require(igraph)){install.packages('igraph'); library(igraph)}
if(!require(sna)){install.packages('sna'); library(sna)}
if(!require(fields)){install.packages('fields'); library(fields)}
if(!require(colorRamps)){install.packages('colorRamps'); library(colorRamps)}
if(!require(locfit)){install.packages('locfit'); library(locfit)}
if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)}
if(!require(gridExtra)){install.packages('gridExtra'); library(gridExtra)}
if(!require(plot3D)){install.packages('plot3D'); library(plot3D)}
if(!require(RcppRoll)){install.packages('RcppRoll'); library(RcppRoll)}

# Loading functions
source(paste(getwd(), '/functions.R', sep=''))
source(paste(getwd(), '/models.R', sep=''))
