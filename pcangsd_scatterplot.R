###### HOW TO MAKE A BEAUTIFUL PCA WITH THE HONEY BEE DATA ##########
# this script reads in .cov files produced from PCAngsd and plots PCA graphs
# written by Daniela Zarate, with help from Erin Calfee 
# last modified: November 17 2020 

# install and load necessary packages
install.packages("RcppCNPy")
install.packages("ggfortify")
library(RcppCNPy)
library(ggfortify)
library(ggplot2)
library(dplyr)

## Read in the .cov files from Documents (the estimated covariance matrix)
honeybees.cov.matrix <- read.table("~/Desktop/ChapterOne-GenomicAdmixture/pcangsd/all.bees.pcangsd.cov", stringsAsFactors = F) 

## Import the Identification List for Legend Colors 
pcangsd.names <- read.table("~/Desktop/ChapterOne-GenomicAdmixture/pcangsd/pcangsd.names", sep = "\n", header = F)

# Use Eigen function on the imported .cov matric 
pca_eigen <- eigen(honeybees.cov.matrix)

# extract eigenvectors = pc's
pcs_eigen <- pca_eigen$vectors
# give the columns names from PC1 to PCX (the end of the columns)
colnames(pcs_eigen) <- paste0("PC", 1:ncol(pcs_eigen))

# normalize eigenvalues to % variance explained
# note: prcomp returns the standard deviations, so they need to be squared
perc_var_pcs_eigen <- 100*pca_eigen$values/sum(pca_eigen$values)

# bind the pcs_eigen matrix to the list matrix to consolidate data 
test <- cbind(pcs_eigen, pcangsd.names)

## create a high resolution figure for export
png("honeybees.pcangsd.jpg", width = 12, height = 7, res = 600, units = 'in', 
    type = 'cairo') # for high resolution figure 

# make a name vector 

names <- c("A. m. scutellata (A)", "A. m. mellifera (M)", "A. m. iberiensis (M)", "A. m. ligustica (C)", "A. m. carnica (C)", "A. m. syriaca (O)", "A. m. anatoliaca (O)", "San Diego (AHB)", "Mexico (AHB)", "Costa Rica (AHB)", "Panama (AHB)")

colors <- c("red", "black", "grey", "yellow", "orange", "aquamarine", "cyan", "orchid", "deeppink",  "tomato", "hotpink4")

legend_title <- "Subspecies / Sample Site"

bee.boo <- pcangsd.names$V1

yg <- ggplot(test, aes(x=PC1, y=PC2, color=bee.boo)) +
  geom_point() +
  scale_color_manual(legend_title, breaks = names, values = colors) + 
  xlab(paste0("PC1 ", round(perc_var_pcs_eigen[1], 2), "%")) +
  ylab(paste0("PC2 ", round(perc_var_pcs_eigen[2], 2), "%"))

# remove ugly grey background, leave grid, make legend labels/text bigger
yg <- yg + theme_bw() + theme(legend.text=element_text(size=14))
# make axes labels bigger 
yg <- yg + theme(axis.title=element_text(size=20))
# make legend title bigger
yg <- yg + theme(legend.title=element_text(size=16))
# make color dots in legend bigger
yg <- yg + guides(colour = guide_legend(override.aes = list(size=8)))

yg
dev.off()

## create a high resolution figure for export
png("honeybees.pcangsd.jpg", width = 12, height = 7, res = 600, units = 'in', 
    type = 'cairo') # for high resolution figure 
