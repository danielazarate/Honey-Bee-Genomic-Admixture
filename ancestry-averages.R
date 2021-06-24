# calculate averages and ranges of admixture results 
# for all populations (all.bees.K4.qopt)

# set working space 
setwd("~/Desktop/ChapterOne-GenomicAdmixture/")

# import the list 
bees.K4 <- read.delim("~/Desktop/ChapterOne-GenomicAdmixture/genomic-admixture-analysis/all.bees.K4.qopt", header = F, sep = " ") 
# remove the random trailing column 
bees.K4 <- bees.K4[,1:4]

# subset the different populations 
African <- bees.K4[1:10,]
west.Euro <- bees.K4[11:30,]
east.Euro <- bees.K4[31:49,]
Mid.East <- bees.K4[50:69,]
San.Diego <- bees.K4[70:84,]
Mexico <- bees.K4[85:99,]
Costa.Rica <- bees.K4[100:114,]
Panama <- bees.K4[115:129,]

# Calucluate Means across columns 
colMeans(Panama[sapply(Panama, is.numeric)])

# calculate standard error of the mean by calculating standard deviation 
# and then dividing by the square root of the sample size. 
sapply(Panama,function(x)sd(x)/sqrt(15))

# calculate min and max to each column of a data frame 
apply(Panama, 2, FUN=min)
