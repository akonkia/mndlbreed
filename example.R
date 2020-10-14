source("breedDragons/makeDragons.R")
library(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#User input
offspring <- 100 # number of offspring to be generated
t <- sample(alleles, 2) #  Number of traits to be taken. Up to 7 traits available.  For now, we assume that one trait goes on one individual chromosome.

#alleles[alleles %in% t]



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
freq <- getFreq(sward)


min <- freq[freq$Count == min(freq$Count),]

paste0("The rarest dragon genotype(s):\n ", paste(min$Genotype, collapse = ", "),"\n was/were spotted ", min$Count[1], " time(s).")
max <- freq[freq$Count == max(freq$Count),]

paste0("In the population of dragon babies, the most often occurring genotype(s):\n ", paste(max$Genotype, collapse = ", ")," was/were spotted ", max$Count[1], " time(s).")



######Get phenotype of a dragon

getPhenotype(sward[[3]])

#Print photypic table
phenoTable

#unique(phenoTable$Phenotype)
collapsed <- collapsePhenoTable(phenoTable)

shortenPhenoTable(t)


type <- c(unlist(strsplit(collapseGenotype(t), " ")))
type <- c(type, tolower(type))
short <- phenoTable$Phenotype[phenoTable$Genotype %in% type]
shortTable <- phenoTable[phenoTable$Phenotype %in% short,]
shortTable <- collapsePhenoTable(shortTable)

# Show frequencies of the phenotypes
showPhenFreq(sward)
