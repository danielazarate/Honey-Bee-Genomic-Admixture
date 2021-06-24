# script to calculate pi for each population 

#note that the "edit" files are pestPG files sans chr 1.0 because I needed to manually remove this chr when i did the jackknives 
# import the file
san <- read.delim("~/Desktop/ChapterOne-GenomicAdmixture/thetastat-pi-calcs/san.diego.pestPG", header = TRUE)
mex <- read.delim("~/Desktop/ChapterOne-GenomicAdmixture/thetastat-pi-calcs/mexico.stat.pestPG", header = TRUE)
cos <- read.delim("~/Desktop/ChapterOne-GenomicAdmixture/thetastat-pi-calcs/costa.rica.stat.pestPG", header = TRUE)
pan <- read.delim("~/Desktop/ChapterOne-GenomicAdmixture/thetastat-pi-calcs/panama.stat.pestPG", header = TRUE)
a.pop <- read.delim("~/Desktop/ChapterOne-GenomicAdmixture/thetastat-pi-calcs/a.pop.pestPG", header = TRUE)
m.pop <- read.delim("~/Desktop/ChapterOne-GenomicAdmixture/thetastat-pi-calcs/m.pop.pestPG", header = TRUE)
c.pop <- read.delim("~/Desktop/ChapterOne-GenomicAdmixture/thetastat-pi-calcs/c.pop.pestPG", header = TRUE)
o.pop <- read.delim("~/Desktop/ChapterOne-GenomicAdmixture/thetastat-pi-calcs/o.pop.pestPG", header = TRUE)


# remove all the rows that are unassigned groups 
# remove the rows with this particular string attachced

# sample populations
san.UnClean <- san[!grepl("GroupUn", san$Chr),]
mex.UnClean <- mex[!grepl("GroupUn", mex$Chr),]
cos.UnClean <- cos[!grepl("GroupUn", cos$Chr),]
pan.UnClean <- pan[!grepl("GroupUn", pan$Chr),]

# reference populations 
a.pop.UnClean <- a.pop[!grepl("GroupUn", a.pop$Chr),]
m.pop.UnClean <- m.pop[!grepl("GroupUn", m.pop$Chr),]
c.pop.UnClean <- c.pop[!grepl("GroupUn", c.pop$Chr),]
o.pop.UnClean <- o.pop[!grepl("GroupUn", o.pop$Chr),]

# San Diego 
# sum the pairwise estimate column and the #nsites column 
san.tP.sum <- sum(san.UnClean$tP)
san.sites.sum <- sum(san.UnClean$nSites)

# Mexico 
mex.tP.sum <- sum(mex.UnClean$tP)
mex.sites.sum <- sum(mex.UnClean$nSites)

# Costa Rica 
cos.tP.sum <- sum(cos.UnClean$tP)
cos.sites.sum <- sum(cos.UnClean$nSites)

#Panama
pan.tP.sum <- sum(pan.UnClean$tP)
pan.sites.sum <- sum(pan.UnClean$nSites)

#Clade A 
a.tP.sum <- sum(a.pop.UnClean$tP)
a.sites.sum <- sum(a.pop.UnClean$nSites)

#Clade M
m.tP.sum <- sum(m.pop.UnClean$tP)
m.sites.sum <- sum(m.pop.UnClean$nSites)

#Clade C
c.tP.sum <- sum(c.pop.UnClean$tP)
c.sites.sum <- sum(c.pop.UnClean$nSites)

#Clade C
o.tP.sum <- sum(o.pop.UnClean$tP)
o.sites.sum <- sum(o.pop.UnClean$nSites)

# calculate pi by dividing the sum of pairwise estimate and the sum of #nsites column 
san.pi <- (san.tP.sum/san.sites.sum)
san.pi = 0.01019705

mex.pi <- (mex.tP.sum/mex.sites.sum)
mex.pi = 0.01158389

cos.pi <- (cos.tP.sum/cos.sites.sum)
cos.pi = 0.01116341

pan.pi <- (pan.tP.sum/pan.sites.sum)
pan.pi = 0.01136556

a.pi <- (a.tP.sum/a.sites.sum)
a.pi = 0.008417279 

m.pi <- (m.tP.sum/m.sites.sum)
m.pi = 0.004383432

c.pi <- (c.tP.sum/c.sites.sum)
c.pi = 0.003463934

o.pi <- (o.tP.sum/o.sites.sum)
o.pi = 0.006428551

# sum the Watterson estimate column

san.tW.sum <- sum(san.UnClean$tW)

mex.tW.sum <- sum(mex.UnClean$tW)

cos.tW.sum <- sum(cos.UnClean$tW)

pan.tW.sum <- sum(pan.UnClean$tW)

a.tW.sum <- sum(a.pop.UnClean$tW)

m.tW.sum <- sum(m.pop.UnClean$tW)

c.tW.sum <- sum(c.pop.UnClean$tW)

o.tW.sum <- sum(o.pop.UnClean$tW)

# calculate pi by dividing the sum of Watterson estimate and the sum of #nsites column 

san.watterson.pi <- (san.tW.sum/san.sites.sum)
san.watterson.pi = 0.01147647

mex.watterson.pi <- (mex.tW.sum/mex.sites.sum)
mex.watterson.pi = 0.01274064

cos.watterson.pi <- (cos.tW.sum/cos.sites.sum)
cos.watterson.pi = 0.01128495

pan.watterson.pi <- (pan.tW.sum/pan.sites.sum)
pan.watterson.pi = 0.01173349

a.watterson.pi <- (a.tW.sum/a.sites.sum)
a.watterson.pi = 0.01111024

m.watterson.pi <- (m.tW.sum/m.sites.sum)
m.watterson.pi = 0.004689856

c.watterson.pi <- (c.tW.sum/c.sites.sum)
c.watterson.pi = 0.003685314

o.watterson.pi <- (o.tW.sum/o.sites.sum)
o.watterson.pi = 0.007138273

# Tajimas's D mean 
san.tajima.sum <- sum(san$Tajima)
san.tajima.pi <- (san.tajima.sum/san.sites.sum)
san.tajima.pi
-1.303429e-05

mex.tajima.sum <- sum(mex$Tajima)
mex.tajima.pi <- (mex.tajima.sum/mex.sites.sum)
mex.tajima.pi
-8.372489e-07

cos.tajima.sum <- sum(cos.UnClean$Tajima)
cos.tajima.pi <- (cos.UnClean.tajima.sum/cos.UnClean.sites.sum)
cos.tajima.pi
-3.378026e-07

pan.tajima.sum <- sum(pan$Tajima)
pan.tajima.pi <- (pan.tajima.sum/pan.sites.sum)
pan.tajima.pi
-4.594793e-07

a.tajima.sum <- sum(a.pop.UnClean$Tajima)
a.tajima.pi <- (a.tajima.sum/a.sites.sum)
a.tajima.pi
-1.936199e-06

m.tajima.sum <- sum(m.pop.UnClean$Tajima)
m.tajima.pi <- (m.tajima.sum/m.sites.sum)
m.tajima.pi

c.tajima.sum <- sum(c.pop.UnClean$Tajima)
c.tajima.pi <- (c.tajima.sum/c.sites.sum)
c.tajima.pi

o.tajima.sum <- sum(o.pop.UnClean$Tajima)
o.tajima.pi <- (o.tajima.sum/o.sites.sum)
o.tajima.pi


### jackknifing it 

############ san.diego ############ 

# Original Pi Calculation using all chromosomes 
san.diego.tP.pi <- 0.00346 
san.diego.tW.pi <- 0.00369 
san.diego.tajima <- -1.303429e-05 

# Group 1 is already excluded upon importing, just remove the Unassigned Groups now 
san.diego.Group1 <- read.delim("~/Documents/SCRIPTS/pi_calculations/san.diego.edit.pestPG", header = TRUE)
# remove unassigned groups
san.diego.Group1 <- san.diego.Group1[!grepl("GroupUn", san.diego.Group1$Chr),]

#san.diego.Group1 <- san.diego.UnClean[!grepl("Group1.", san.diego.UnClean$Chr),]
san.diego.Group1.Tajima.sum <- sum(san.diego.Group1$Tajima)
san.diego.Group1.sites.sum <- sum(san.diego.Group1$nSites)
san.diego.Group1.pi <- (san.diego.Group1.Tajima.sum/san.diego.Group1.sites.sum)
san.diego.Group1.pi


# group 2 
san.diego.Group2 <- san.diego.Group1[!grepl("Group2", san.diego.Group1$Chr),]
san.diego.Group2.Tajima.sum <- sum(san.diego.Group2$Tajima)
san.diego.Group2.sites.sum <- sum(san.diego.Group2$nSites)
san.diego.Group2.pi <- (san.diego.Group2.Tajima.sum/san.diego.Group2.sites.sum)
san.diego.Group2.pi

# group 3 
san.diego.Group3 <- san.diego.Group2[!grepl("Group3", san.diego.Group2$Chr),]
san.diego.Group3.Tajima.sum <- sum(san.diego.Group3$Tajima)
san.diego.Group3.sites.sum <- sum(san.diego.Group3$nSites)
san.diego.Group3.pi <- (san.diego.Group3.Tajima.sum/san.diego.Group3.sites.sum)
san.diego.Group3.pi 


# group 4 
san.diego.Group4 <- san.diego.Group3[!grepl("Group4", san.diego.Group3$Chr),]
san.diego.Group4.Tajima.sum <- sum(san.diego.Group4$Tajima)
san.diego.Group4.sites.sum <- sum(san.diego.Group4$nSites)
san.diego.Group4.pi <- (san.diego.Group4.Tajima.sum/san.diego.Group4.sites.sum)
san.diego.Group4.pi

# group 5 
san.diego.Group5 <- san.diego.Group4[!grepl("Group5", san.diego.Group4$Chr),]
san.diego.Group5.Tajima.sum <- sum(san.diego.Group5$Tajima)
san.diego.Group5.sites.sum <- sum(san.diego.Group5$nSites)
san.diego.Group5.pi <- (san.diego.Group5.Tajima.sum/san.diego.Group5.sites.sum)
san.diego.Group5.pi

# group 6 
san.diego.Group6 <- san.diego.Group5[!grepl("Group6", san.diego.Group5$Chr),]
san.diego.Group6.Tajima.sum <- sum(san.diego.Group6$Tajima)
san.diego.Group6.sites.sum <- sum(san.diego.Group6$nSites)
san.diego.Group6.pi <- (san.diego.Group6.Tajima.sum/san.diego.Group6.sites.sum)
san.diego.Group6.pi


# group 7 
san.diego.Group7 <- san.diego.Group6[!grepl("Group7", san.diego.Group6$Chr),]
san.diego.Group7.Tajima.sum <- sum(san.diego.Group7$Tajima)
san.diego.Group7.sites.sum <- sum(san.diego.Group7$nSites)
san.diego.Group7.pi <- (san.diego.Group7.Tajima.sum/san.diego.Group7.sites.sum)
san.diego.Group7.pi

# group 8 
san.diego.Group8 <- san.diego.Group7[!grepl("Group8", san.diego.Group7$Chr),]
san.diego.Group8.Tajima.sum <- sum(san.diego.Group8$Tajima)
san.diego.Group8.sites.sum <- sum(san.diego.Group8$nSites)
san.diego.Group8.pi <- (san.diego.Group8.Tajima.sum/san.diego.Group8.sites.sum)
san.diego.Group8.pi

# group 9 
san.diego.Group9 <- san.diego.Group8[!grepl("Group9", san.diego.Group8$Chr),]
san.diego.Group9.Tajima.sum <- sum(san.diego.Group9$Tajima)
san.diego.Group9.sites.sum <- sum(san.diego.Group9$nSites)
san.diego.Group9.pi <- (san.diego.Group9.Tajima.sum/san.diego.Group9.sites.sum)
san.diego.Group9.pi

# group 10 
san.diego.Group10 <- san.diego.Group9[!grepl("Group10", san.diego.Group9$Chr),]
san.diego.Group10.Tajima.sum <- sum(san.diego.Group10$Tajima)
san.diego.Group10.sites.sum <- sum(san.diego.Group10$nSites)
san.diego.Group10.pi <- (san.diego.Group10.Tajima.sum/san.diego.Group10.sites.sum)
san.diego.Group10.pi

# group 11 
san.diego.Group11 <- san.diego.Group10[!grepl("Group11", san.diego.Group10$Chr),]
san.diego.Group11.Tajima.sum <- sum(san.diego.Group11$Tajima)
san.diego.Group11.sites.sum <- sum(san.diego.Group11$nSites)
san.diego.Group11.pi <- (san.diego.Group11.Tajima.sum/san.diego.Group11.sites.sum)
san.diego.Group11.pi

# group 12 
san.diego.Group12 <- san.diego.Group11[!grepl("Group12", san.diego.Group11$Chr),]
san.diego.Group12.Tajima.sum <- sum(san.diego.Group12$Tajima)
san.diego.Group12.sites.sum <- sum(san.diego.Group12$nSites)
san.diego.Group12.pi <- (san.diego.Group12.Tajima.sum/san.diego.Group12.sites.sum)
san.diego.Group12.pi

# group 13 
san.diego.Group13 <- san.diego.Group12[!grepl("Group13", san.diego.Group12$Chr),]
san.diego.Group13.Tajima.sum <- sum(san.diego.Group13$Tajima)
san.diego.Group13.sites.sum <- sum(san.diego.Group13$nSites)
san.diego.Group13.pi <- (san.diego.Group13.Tajima.sum/san.diego.Group13.sites.sum)
san.diego.Group13.pi

# group 14
san.diego.Group14 <- san.diego.Group13[!grepl("Group14", san.diego.Group13$Chr),]
san.diego.Group14.Tajima.sum <- sum(san.diego.Group14$Tajima)
san.diego.Group14.sites.sum <- sum(san.diego.Group14$nSites)
san.diego.Group14.pi <- (san.diego.Group14.Tajima.sum/san.diego.Group14.sites.sum)
san.diego.Group14.pi

# group 15
san.diego.Group15 <- san.diego.Group14[!grepl("Group15", san.diego.Group14$Chr),]
san.diego.Group15.Tajima.sum <- sum(san.diego.Group15$Tajima)
san.diego.Group15.sites.sum <- sum(san.diego.Group15$nSites)
san.diego.Group15.pi <- (san.diego.Group15.Tajima.sum/san.diego.Group15.sites.sum)
san.diego.Group15.pi

# vector of all jackknife mean estimates
san.diego.Tajima.jackknives <- c(san.diego.tajima, san.diego.Group1.pi, san.diego.Group2.pi, san.diego.Group3.pi, san.diego.Group4.pi, san.diego.Group5.pi, san.diego.Group6.pi, san.diego.Group7.pi, san.diego.Group8.pi, san.diego.Group9.pi, san.diego.Group10.pi, san.diego.Group11.pi, san.diego.Group12.pi, san.diego.Group13.pi, san.diego.Group14.pi, san.diego.Group15.pi)

# mean of all the means of all jackknives 
mean(san.diego.Tajima.jackknives)
-1.655431e-06
# sum of squares # obtain from online calculator 
# got to find an R equivalent ...
sumofsq <- 1.3826268286945E-10
# (n-1)/16 
n <- (16 - 1) / 16
# multiple the sum of squares by (n-1)/16
oi <- n*sumofsq
# square root 
sqrt(oi)
## standard error 
0.000543
