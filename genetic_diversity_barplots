#!/bin/env/R 
# Script for plotting barplot of genetic diversity measurements in honey bee pops. 
# Author: Daniela Zarate, PhD.
# Built and run under R.4.3.3 - Angel Food Cake
#_______________________________________________________________________________#

# Load some packages
library(ggplot2)
library(dplyr)

#_______________________________________________________________________________#

# Read in the data 
diversity <- read.csv("~/Desktop/HONEY_BEE_ADMIXTURE/diversity.measures.csv")
attach(diversity)

# Quick barplot for rough visual:
barplot(diversity$Pairwise.Estimator)

# plots site alphabetically, to plot temporally, specify factor levels: 
diversity$Population  <-  factor(
  diversity$Population,
  levels=c("African", "West.Europe", "East.Europe", "Middle.Eastern",
           "San.Diego", "Mexico", "Costa.Rica", "Panama"))
# Check the levels to make sure they are in the correct order:
levels(diversity$Population)

# Make a vector of all the standard errors per pop:
testvec <- c(0.0018, 0.0007, 0.0003, 0.0005, 0.0018, 0.0023, 0.0024, 0.0024)

# Add the vector as a new column to the diversity df:
diversity$sem <- testvec

# Basic barplot of Watterson Estimator: 
ggplot(data=diversity, aes(x=Population, y=Pairwise.Estimator, fill=Population)) +
  geom_bar(stat="identity") + theme_bw() +
  geom_errorbar(aes(x=Population, ymin=Pairwise.Estimator-sem, ymax=Pairwise.Estimator+sem))+
  scale_fill_manual(values=c("red", "black", "yellow", "cyan","orchid", "deeppink", "tomato","hotpink4"),
      name = "Population", labels = c("African", "Western Europe", 
   "Eastern Europe", "Middle Eastern", "San Diego", "Mexico", "Costa Rica", "Panama")) +
  ylab("Genome-wide pairwise estimator") + xlab(" ") + theme(axis.title.y = element_text(size = 16)) +
 scale_x_discrete(labels=c("", "", "", "","", "", "", "" )) + theme(axis.ticks.x=element_blank())

#_______________________________________________________________________________#

# Set the high-resolution png file:
png(("~/Genome_Pairwise_Estimator.png"),
    width = 10, height = 7, res = 300,units = 'in', )
dev.off()
