source("breedDragons/makeDragons.R")
library(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#User input
offspring <- 100 # number of offspring to be generated
t <- sample(alleles, 7) #  Number of traits to be taken. Up to 7 traits available.  For now, we assume that one trait goes on one individual chromosome.

# Start with creating two dragon parents
drag1 <- createDragon(t)
drag2 <- createDragon(t)

# Plot genotypes of the parents
showPair(drag1, drag2)

# create a sward of dragons, offspring (n) of the selected parents. Outputs plot under 4 traits, and table over 3 traits. Need to modify.
sward <- makeSward(drag1, drag2, offspring)

# Create summary of the sward genotypes
showGenFreq(sward)

#Make a plot of the sward genotypes
plotGenFreq(sward)

#Create a table of the sward genotypes
dim(getFreq(sward))

# Plot a selected dragon from the sward. You can add a name as a second argument.
plotGenotype(sward[[2]], "Elvis")

# Show frequencies of the genotypes (table)
freq <- getFreq(sward)

# Get phenotype of a dragon

getPhenotype(sward[[3]])

# Print the phenotype table
phenoTable

# Collapse phenotype table (ex. put Aa and AA in the same cell, related to the same phenotype)
collapsePhenoTable(phenoTable)

# Create a short collapsed phenotype table, containing only traits chosen for the simulation
shortTable <- shortenPhenoTable(t)

# Show frequencies of the phenotypes
showPhenFreq(sward)
