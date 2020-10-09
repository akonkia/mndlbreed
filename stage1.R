source("makeDragons.R")
library(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#User input
offspring <- 100 # number of offspring to be generated. Slider with options 1-10 000, for ex.
t <- sample(alleles, 2) #  Number of traits to be taken. Up to 7 traits available.  For now, we assume that one trait goes on one individual chromosome.
# Add option to pick alleles for parents. Ex. radio buttons for every option. Would need to be interactively generated as step 2. For now it is random

# Start with creating two dragon parents
drag1 <- createDragon(t)
drag2 <- createDragon(t)

# Plot genotypes of the parents
showParents(drag1, drag2)

#create a sward of dragons, offspring (n) of the selected parents
sward <- makeSward(drag1, drag2, offspring)

# Plot a selected dragon from the sward
plotGenotype(sward[[2]])

# Show frequencies of the genotypes
showGenFreq(sward)

######Get phenotype of a dragon
getPhenotype(sward[[3]])

#Print photypic table
phenoTable

# Show frequencies of the phenotypes
showPhenFreq(sward)
