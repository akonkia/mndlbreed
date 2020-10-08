source("makeDragons.R")
library(ggpubr)
#library(RColorBrewer)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#User input
n <- 2 # number of offspring to be generated. Slider with options 1-10 000, for ex.
t <- sample(alleles, 3) #  Number of traits to be taken. Up to 7 traits available.  For now, we assume that one trait goes on one individual chromosome.
# Add option to pick alleles for parents. Ex. radio buttons for every option. Would need to be interactively generated as step 2. For now it is random

drag1 <- createDragon(t)
drag2 <- createDragon(t)
p1 <- plotGenotype(drag1, "Dragon 1")

getGenotype(drag1)
p2 <- plotGenotype(drag2, "Dragon 2")
getGenotype(drag2)

showParents(drag1, drag2)

collapseGenotype(drag2)

sward <- list()
for (i in 1:n){
  sward <- list.append(sward, breedDragon(drag1, drag2))
}

getGenotype(sward[1]) #vector with ordered genotype
plotGenotype(sward[[2]])

# title <- data.frame(cbind(1:3, 1:6))
# title1 <- ggplot(data = title)+
#   theme_void()+
#   annotate(geom = "text", x = 1, y = max(as.numeric(df$Length))+0.5, label = c("Parent 1"))
# ggarrange(title1, p2)









##TODO adapt following plot to show genotypes of parents, and then offspring


#plot genes on chromosomes
chrom_key2 <- setNames(object = as.character(1:5), 
                      nm = df$chroms)
df <- subset(chrom_sizes, chroms %in% genes_bands$chromosome)

rownames(df) <- 1:5
#windows(20,10)
ggplot(data = df)+
  theme_void(base_size = 12)+
  geom_rect(aes(xmin = as.numeric(rownames(df)) - 0.2, 
                xmax = as.numeric(rownames(df)) + 0.2, 
                ymax = Length, 
                ymin = 0), 
                fill = "light blue")+
  # rotate the plot 90 degrees
  coord_flip()+
  scale_x_discrete(name = "Chromosome", limits = df$chroms)+
# geom_rect(data = genes_bands, aes(xmin = as.numeric(chromosome) - 0.2, 
#                                   xmax = as.numeric(chromosome) + 0.2, 
#                                   ymax = end, ymin = start, fill = rev(genes_bands$genes)))+ #remove label, sort order, make that a gene
  geom_text(data = genes, #add subsetting for genes that were chosen, ex data = subset(sample_cns, sample_cns$CNA == "gain")
                  aes(x = chromosome, y = start, label = Gene), 
                  color = c("red", "blue", "black", "purple", "dark green"), show.legend = FALSE,
                  size = 6)





###################Dragon chromosomes###########################################
# Create chromosome database
#chromosomes = c("chr1", "chr2", "chr3", "chr4","chr5", "chr6", "chr7")
#chroms = c("chr1", "chr2", "chr3", "chr4","chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12")



#chrom_sizes <- data.frame(chroms, length)

#chrom_sizes <- arrange(chrom_sizes, desc(length))

#chrom_sizes$chroms <- factor(chrom_sizes$chroms,  levels = chroms)

###############################################################################

# chrom_key <- setNames(object = as.character(c(1, 2, 3, 4, 5, 6, 7)), 
#                       nm = chroms)




# 
# #genes$Chromosome <- factor(genes$Chromosome,  levels = chromosome)
# 
# ggplot(df, aes(y = Length, x = Chromosome,  fill = Set))+
#   theme_minimal(base_size = 12)+
#   geom_col(position = "dodge2")+
#   theme(legend.position = "none")+
#   geom_text(aes(y = Start, x = Chromosome, label = Chromosome))
# 
