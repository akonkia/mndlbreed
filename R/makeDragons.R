library("scales") # for axis labels notation#
#library("dplyr")
library(rlist)
library(tidyverse)


###################################Dragon genes################################
gene_names <- c("Body colour", "Horn", "Wings", "Tail colour", "Tail spikes", "Toes number", "Fire breathing")

#alleles = data.frame(c("A", "a"), c("B", "b"), c("C", "c"), c("D", "d"),c("E", "e"))
alleles <- data.frame(c("A", "a"), c("H", "h"), c("W", "w"), c("T", "t"),c("S", "s"), c("M", "m"),c("F", "f"))
colnames(alleles) <- c("Body and head color", "Horn", "Wing color", "Tail color", "Number of tail spikes", "Number of toes", "Fire-breathing")
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

################################Phenotyping table#############################
Genotype <- c("FF", "Ff", "ff", "MM", "Mm", "mm", "SS", "Ss", "ss", "TT", "Tt", "tt", "WW", "Ww", "ww", "HH", "Hh", "hh", "AA", "Aa", "aa")

Phenotype <- c("Fire-breathing", "Fire-breathing", "Does not breath fire", "Four toes", "Four toes", "Three toes", "Five spikes on tail", "Five spikes on tail", "Four spikes on tail",
               "Red tail", "Red tail", "Yellow tail", "Red wings", "Red wings", "Yellow wings", "Horn", "Horn", "No horn", "Blue body and head", "Blue body and head", "Green body and head")

phenoTable <- data.frame(Genotype, Phenotype)

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
createDragon <- function(traits){
  n <- 2 #diploidy
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


# ##############Set parents###############
# setParents <- function(drag1, drag2){
#   parents <- list(drag1, drag2)
#   return(parents)
# }
# ########################################

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
  parents <- list(first, second)
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
#helper function to getGenotype
sortGenotype <- function(x) {
  tmp <- sort(x)
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
  sorted <- sortGenotype(gen)
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
    
    geom_label(data = df, 
               aes(x = Chromosome, y = Start, label = allele), 
               show.legend = FALSE,
               fill = colors[1:dim(df)[1]], 
               colour = "white",
               fontface = "bold",
               family = "mono",
               size = 9,
               
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
showPair <- function(dragon, dragon2){
  
  p1 <- plotGenotype(dragon, "Dragon 1")
  p2 <- plotGenotype(dragon2, "Dragon 2")
  p3 <- grid.arrange(p1, p2, ncol = 1)
  return (p3)
}


###############Collapse genotype############
###Sort alleles, so that aA or Aa is presented as Aa
collapseGenotype <- function(dragon) {
  gen <- getGenotype(dragon)
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

##################Output frequency plot/table for genotype########

plotGenFreq <- function(sward){
  
  freq <- getFreq(sward)
  if (dim(freq)[1] > 100){
    freq <- freq[freq$Count>1,]
  }
  
  
  colors <- c("#8DA0CB", "#FC8D62", "#66C2A5", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
  colors <- rep(colors,200)
  output <- ggplot(freq, aes(x = Genotype, y = Count)) +
    theme_minimal()+
    
    geom_bar(fill = colors[1:dim(freq)[1]], stat = "identity") +
    geom_text(aes(label = Count), hjust = -0.5)+
    
    theme(panel.grid.major = element_blank(), axis.text.x = element_blank(), panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.title.y=element_blank(),axis.title.x=element_blank())+
    coord_flip()+
    ggtitle ("How many dragons had the same genotype?")

  return(output)
}

###################Count freq of the genotypes#############
countFreq <- function(sward){
  for (i in 1:length(sward)){
    if (i == 1){
      pool <- data.frame(collapseGenotype(sward[[i]]), t(getGenotype(sward[[i]])))
      colnames(pool) <- c("Genotype", 1:(length(getGenotype(sward[[i]]))))
    }else{
      pool <- rbind(pool, c(collapseGenotype(sward[[i]]), getGenotype(sward[[i]])))
    }
  }
  return(pool)
}

#trick to turn list to data.frame, strings as factors:
#df <- data.frame(matrix(unlist(list), nrow=length(list), byrow=T))
#strings not as factors:
#df <- data.frame(matrix(unlist(l), nrow=132, byrow=T),stringsAsFactors=FALSE)

####################create sward of dragons#############
makeSward <- function(drag1, drag2, offspring){
  sward <- list()
  for (i in 1:offspring){
    sward <- list.append(sward, breedDragon(drag1, drag2))
  }
  return(sward)
}


######create phenotype###################################
getPhenotype <- function(dragon){
  type <- c(unlist(strsplit(collapseGenotype(dragon), " ")))
  phenotype <- phenoTable[phenoTable$Genotype %in% type,]
  rownames(phenotype) <- NULL
  return(phenotype)
}

########Phenotype sward############################
#####Phenotypes every dragon in a sward###########
phenotypeSward <- function(sward){
  phenotypes <- list()
  for (each in sward){
    type <- c(unlist(strsplit(collapseGenotype(each), " ")))
    tmp <- paste(phenoTable[phenoTable$Genotype %in% type,]$Phenotype, collapse = ', ')
    phenotypes <- list.append(phenotypes, tmp)
  }
  df <- data.frame(matrix(unlist(phenotypes), nrow=length(phenotypes), byrow=T), rep(1, length(phenotypes)))
  colnames(df) <- c("Phenotype", "Count")
  return(df)
}

############## Show frequencies of the phenotypes
showPhenFreq <- function(sward){
  phenotypes <- phenotypeSward(sward)
  freq <- count (phenotypes, phenotypes$Phenotype)
  freq <- freq[order(-freq$n),]
  colnames(freq) <- c("Phenotype", "Count")
  rownames(freq) <- NULL
  return(freq)
}


##############Get genotype freq table
getFreq <- function(sward){
  pool <- countFreq(sward)
  freq <- count (pool, pool$Genotype)
  #sort the freq table
  freq <- freq[order(nrow(freq):1),]
  rownames(freq) <- 1:dim(freq)[1]
  colnames(freq) <- c("Genotype", "Count")
  freq$Genotype <- factor(freq$Genotype)
  return(freq)
}

#################Give plot or table of genfreq###
showGenFreq <- function(sward){
  
  freq <- getFreq(sward)
  
  if (dim(freq)[1] <9){
    
    output <- plotGenFreq(sward)
  }else{
    output <- freq
  }
  
  return (output)
}

##########Collapse phenoTable ####
collapsePhenoTable <- function(phenoTable){
  
  temp <- phenoTable[!duplicated(phenoTable$Phenotype),]
  dup <- phenoTable[duplicated(phenoTable$Phenotype),]
  
  j= 1
  for (i in 1:dim(temp)[1]){
    if(temp$Phenotype[i] %in% dup$Phenotype){
      temp$Genotype[i] <- paste(temp$Genotype[i], dup$Genotype[j])
      j <- j + 1
    }
  }
  rownames(temp) <- NULL
  return (temp)
}


############Shorten phenoTable based on traits chosen######
shortenPhenoTable <- function(traits){
  type <- c(unlist(strsplit(collapseGenotype(traits), " ")))
  type <- c(type, tolower(type))
  short <- phenoTable$Phenotype[phenoTable$Genotype %in% type]
  shortTable <- phenoTable[phenoTable$Phenotype %in% short,]
  shortTable <- collapsePhenoTable(shortTable)
  rownames(shortTable) <- NULL
  
  return(shortTable)
}


homozygoteRe <- function(dragon){
  dragon$set1 <- tolower(dragon$set1)
  dragon$set2 <- tolower(dragon$set2)
  return(dragon)
}

homozygoteDo <- function(dragon){
  dragon$set1 <- toupper(dragon$set1)
  dragon$set2 <- toupper(dragon$set2)
  return(dragon)
}

heterozygote <- function(dragon){
  dragon$set1 <- toupper(dragon$set1)
  dragon$set2 <- tolower(dragon$set2)
  return(dragon)
}
