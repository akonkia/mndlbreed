library("scales") # for axis labels notation#
#library("dplyr")
library(rlist)
library(tidyverse)


###################################Dragon genes################################
gene_names <- c("Body colour", "Horn", "Wings", "Tail colour", "Tail spikes", "Toes number", "Fire breathing")

#alleles = data.frame(c("A", "a"), c("B", "b"), c("C", "c"), c("D", "d"),c("E", "e"))
alleles <- data.frame(c("A", "a"), c("H", "h"), c("W", "w"), c("T", "t"),c("S", "s"), c("M", "m"),c("F", "f")) 
length = c(100L, 150L, 180L, 200L, 250, 280L, 300L)
chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7")
start = c(50,80, 40, 100,150,80,260)
end = c(60,90, 50,110,160,90,270)

genes <- data.frame(chromosome, length, start, end, gene_names, t(alleles))
rownames(genes) <-  NULL
genes$chromosome <- factor(genes$chromosome,  levels = chromosome)
colnames(genes) <- c("Chromosome", "Length", "Start", "End", "Gene", "Allele 1", "Allele 2")

# save(genes, file = "genes.RData")
# load("genes.RData")

##############################################################################


##################create haplotype####
createHaplotype <- function(traits){
  pick <- c()
  t <- traits
  for (i in 1:length(t)){
    pick <- c(pick, sample(t[[i]], 1))
  }
  return(pick)
}
######################################

#################create dragon########
createDragon <- function(traits, diploidy = 2){
  n <- diploidy
  t<- traits
  dragon <- list()
  for (i in 1:n){
    set <- createHaplotype(t)
    dragon <- list.append(dragon, set)
  }
  dragon <- data.frame(dragon)
  colnames(dragon)<- c("set1", "set2") #will throw an error for diploidy different than 2
  return(dragon)
}
#######################################


##############Set parents###############
setParents <- function(drag1, drag2){
  parents <- list(drag1, drag2)
  return(parents)
}
########################################

###############Generate gametes#########
generateGamete <- function(dragon){
  traits <- dim(dragon)[1]
  gamete <- vector(length = traits)
  for (i in 1:traits){
    allele <- unlist(sample(dragon[i,],1))
    gamete[i] <- allele
  }
  return (gamete)
}
#######################################

#################Breed dragons#########
breedDragon <- function(first, second){
  parents <- setParents (drag1, drag2)
  n <- 2
  dragon <- list()
  for (i in 1:n){
    gamete <- generateGamete(parents[[i]])
    dragon <- list.append(dragon, gamete)
  }
  dragon <- data.frame(dragon)
  colnames(dragon)<- c("set1", "set2") #will throw an error for diploidy different than 2
  return(dragon)
}

########################################
###############Sort genotype############
###Sort alleles, so that aA or Aa is presented as Aa
sortGenotype <- function(x) {
  tmp <- x
g <- vector(length = length(tmp))

for (char in 1:(length(tmp)-1)){
  
  if ( tmp[char] == tmp[char+1]){
    g[char] <- tmp[char]
    g[char+1] <- tmp[char+1]
  }else if(tolower(tmp[char]) == tolower(tmp[char+1])){
    if (str_detect(tmp[char], "^[:upper:]+$")){
      g[char] <- tmp[char]
      g[char+1] <- tmp[char+1]
    }else{
      g[char] <- tmp[char+1]
      g[char+1] <- tmp[char]
    }
  }
}

return (g)}
########################################

###############Get Genotype############
##Present genotype in a succinct manner
getGenotype <- function(dragon){
  gen <- unlist(dragon)
  sorted <- sortGenotype(sort(gen))
  return(sorted)
}


####################Plot chromosomes#####

plotChroms <- function(df, allele_set){
  colors <- c("#8DA0CB", "#FC8D62", "#66C2A5", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
  if (allele_set == 1){
    allele = df$`Allele 1`
    pal = "light blue"
  }else{
    allele = df$`Allele 2`
    pal = "light green"
  }
  
  ggplot(data = df)+
    theme_void(base_size = 12)+
    
    geom_rect(aes(xmin = as.numeric(rownames(df)) - 0.2, 
                  xmax = as.numeric(rownames(df)) + 0.2, 
                  ymax = Length, ymin = 0), 
              fill = pal)+
    scale_x_discrete(name = "Chromosome", limits = df$Chromosome)+
    # geom_text(data = df, 
    #           aes(x = 2, y = Length, label = "df"),
    #           colour = "black",
    #           size = 9)+
    
    geom_label(data = df, 
              aes(x = Chromosome, y = Start, label = allele), 
               show.legend = FALSE,
              fill = colors[1:dim(df)[1]], 
              colour = "white",
              fontface = "bold",
              family = "mono",
              size = 9,
              
              #color=colors[1:dim(dragon)[1]]
  )
}
#################Plot genotype#############
plotGenotype <- function(dragon, name = "Dragon"){
  geno <- getGenotype(dragon)

  #create a df with chromosomes containing genes taken under account
  df <- genes[genes$`Allele 1`%in% geno | genes$`Allele 2` %in% geno, ] 
  df$`Allele 1` <- dragon$set1
  df$`Allele 2` <- dragon$set2
  rownames(df) <- 1:dim(df)[1]
  #Set 1
  set1 <- plotChroms(df, 1)
  #Set 2
  set2 <- plotChroms(df, 2)
  plot <- grid.arrange(set1, set2, nrow = 1)
  text1 <- text_grob(name, face = "bold", family = "mono", size = 20)
  text2 <- text_grob(collapseGenotype(geno), family = "mono", size = 16)
  text <- grid.arrange(text1, text2, nrow = 2)
  #title <- str_c(name, collapseGenotype(geno), sep = "\n")
  #text <- text_grob(title, family = "mono", size = 20)
  plot2 <- grid.arrange(text, plot, ncol = 1, heights = c(2/10, 8/10))

  return(plot2)
}
#

###############Show parental genotypes####
showParents <- function(dragon, dragon2){
  
  p1 <- plotGenotype(dragon, "Dragon 1")
  p2 <- plotGenotype(dragon2, "Dragon 2")
  p3 <- grid.arrange(p1, p2, ncol = 1)
  return (p3)
}


###############Collapse genotype############
###Sort alleles, so that aA or Aa is presented as Aa
collapseGenotype <- function(x) {
  gen <- getGenotype(x)
  add <- 0
  short <- ""
  for(i in 1:(length(gen)/2)){
    tmp <- str_c(gen[i+add], gen[i+add+1])
    add <- add + 1
    if (i == 1){
      short <- str_c(short, tmp, sep = "")
    }else{
      short <- str_c(short, tmp, sep = " ")
    }
  }
  return (short)
}