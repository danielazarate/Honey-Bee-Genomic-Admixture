# This script plots ancestry proportions for honey bees, using .qopt file produced from NGSADMIX ancestry estimation program. 

# Created: 7.30.2019 
# by daniela zarate 
##############################################################################

# Import the QOPT file and a corresponding sample name list 

# set working directory
setwd("~/Desktop/ChapterOne-GenomicAdmixture/")


# the library is open...because reading is WHAT? Fundamental. 

library(grid)
library(gridExtra)

# paths to input and output directories
# must run this in the directory all files live, otherwise, change to absolute path
qopt_path = paste0("./balance.bees.2.beagle.gz.K6.qopt")
list_path = paste0("./balance.bees.2.list")

# must create the directory in the same folder
results_path = paste0("./results/")

# read in the QOPT file output from NGSADMIX
ancestry <- read.table(qopt_path, 
stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("C", "O", "M", "A", "Z", "X"))

# download list, replace any full paths using search&replace in vscode
# read in the sample names, use the list of bam files provided for ANGSD

samples <- read.table(list_path, stringsAsFactors = F, header = F)

# convert the list of sample names into a vector in order to label the X-axis 

names.vec = samples$V1

# set up a high resolution PNG file, outputting to results directory. 
png(paste0(results_path, "/balanced.bees.K6.jpg"),
    width = 12, 
    height = 7, 
    res = 300, 
    units = 'in', 
    type = 'cairo')

# create barplot

barplot(t(as.matrix(ancestry)),
        col=c("yellow", "red", "black", "cyan", "green", "pink"), 
        cex.names = 0.5, 
        las = 2, # rotate x axis labels vertical 
        ylab = "ancestry proportion", 
        names.arg = names.vec,
        main = "Genomic Admixture in Africanized honey bees across the americas [NGSAdmix, K = 6] ") 

dev.off() # closes all windows 
################################################################################
