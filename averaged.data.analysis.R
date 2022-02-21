#/bin/env/R 
# This script pertains to the data management and analysis of Chapter III in \
# the PhD dissertation of Daniela Zarate at UC San Diego. 
# Data Summary: Here, we analyze a honey bee behavioral dataset collected \
# from May 2021 - November 2021. Several honey bee colonies were assessed \
# across two different sites (BFS: a  managed apiary and ECR: an unmanaged 
# feral site) over several metrics of nest defense behavior. \

# In this script, we visualize the data using both line graphs and boxplots \
# and then run a repeated measures ANOVA to assess differences in nest defense \
# behavior between sites across months. 

#_______________________________________________________________________

## load libaries used in this analysis
## visualizing the data 
library(lattice)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(tidyr)
library(nlme)
library(MANOVA.RM)
library(lme4)
library(lmerTest)
#_______________________________________________________________________

# initial data importation and preliminary checks
# set working directory 
setwd("~/Documents/honeybee_nestdefense/behavioral_data/colony_defense/")

# read in the data 
data <- read.csv("Consolidated_Monthly_Averages.csv")

# attach and review the structure of the data 
attach(data)
str(data)

# chop off last random column 
data <- data[,1:12]
#_______________________________________________________________________

# Data Visualization 

# Cnvert Month to a factor
data$Month <- as.factor(data$Month)

# plots month alphabetically, to plot temporally, specify factor levels 
Time  <-  factor(data$Month, levels=c("MAY_2021", "JULY_2021", "AUGUST_2021", "OCTOBER_2021",
                                      "NOVEMBER_2021"))
# plot raw data 
xyplot(Stings.on.Flag.cont ~ Time, data = data)

# plot with lines, with different color per colony 
xyplot(Stings.on.Flag.cont ~ Time, group = Colony, data = data, type = "b")


## separate by Site, raw data  
xyplot(Stings.on.Flag.cont ~ Time | Site, data = data)


# weirdnesss going on, check Site column 
unique(Site)
# yup, blank spaces after BFS/ECR making them different...

# plot with lines now 
xyplot(Stings.on.Flag.cont ~ Time | Site, group = Colony, data = data, type = "b")

## ggplot 
base.figure <- ggplot(data, aes(x = Time, y = Stings.on.Flag.cont, colour = Colony)) + geom_point()
base.figure # just points, each colony a colour 

## How to make colony a discrete color scale?
lapply(data$Colony, class) # colony is an integer, let's make it a character
data$Colony <- as.factor(data$Colony)
# doesn't do anything

# add line graph on top of base figure (points), each colony a color 
line.figure <- base.figure + geom_line(aes(group = Colony))
# aes(color = Colony, group = Colony) the "Colour = Colony" doesn't do shit. 

line.figure # one figure, each colony a color


# separating by site on one figure 
site.figure <- ggplot(data, aes(x = Time, y = Stings.on.Flag.cont,
                  group = Colony)) + geom_line(aes(color=Site)) 

site.figure <- site.figure+scale_color_manual(values=c("cornflowerblue", "deeppink", "green", "cyan1"))
site.figure

# figure 1, all colonies, seperated by site for stings on flag 
fig1 <- site.figure  + ggtitle("Stings on Flag Across Month and Between Sites") + 
  theme(plot.title = element_text(hjust = 0.5))  + theme_bw() # produce base 
fig1 <- fig1 + theme(plot.title = element_text(hjust = 0.5)) # center title 
fig1 # line graph, one figure, site specified by color

 ## separate above by site, each site one color 
site.sep.figure <- site.figure + facet_wrap( ~ Site, labeller = label_both) + theme_bw()
site.sep.figure
## separate by site, each colony one color 
site.colony.figure <- line.figure + facet_wrap(~Site) + theme_bw() + ggtitle("Stings on Flag Across Month and Between Sites")
site.colony.figure


## display the mean value per group
# we specify here we average over treatment using the mean function
# start with site figure
mean.site.figure <- site.figure + stat_summary(aes(group = Site, color = Site),
                       geom = "line", fun = mean, size = 3)
# final figure 
mean.site.figure + theme_bw() +
  ggtitle("Stings on Flag Across Month and Between Sites with Averages")
# erase the grey background and add title 

## if I want to specify the colors. 
x1 <- x + scale_colour_manual(name = "Group",
                      values = c("BFS" = "cornflowerblue", "mean BFS" = "blue",
                                 "ECR" = "magenta", "mean ECR" = "hotpink"))


## export figures
# Set the high-resolution png file 
png(("figures/Stings.on.Flag.Across.Month.and.Between.Sites.averages.jpg"),
    width = 10, 
    height = 7, 
    res = 300, 
    units = 'in', 
    )

#_______________________________________________________________________

# More data visualization, this time boxplots 
# Plot boxplots of DRV across month and site, and add average overlaid on plot 

# subset data to select only Month, Julian.Date, Colony, Site, and 1 DRV
# DRV = dependent response variable 
data.sub <- data[,c(1,4,5,10)]
# specift by column name 
data.sub <- data[,c("Month", "Site", "Colony", "Stings.on.Flag.cont")]

# transpose rows (months) into columns by melting 
melt.data <- data.sub %>% 
  spread(Month, Stings.on.Flag.cont)
# observations: missing data, 27 colonies total 
# this introduces NAs for missing data 

# transform data into long format in prep for analysis
# month converted to rows "Month" and DRV into "Stings"
# convert id and time into factor variables
long.data <- melt.data %>%
  gather(key = "Month", value = "Stings", MAY_2021, JULY_2021, AUGUST_2021, OCTOBER_2021, NOVEMBER_2021) %>%
  convert_as_factor(Colony, Month)
# check out the data 
head(long.data)


# Compute some summary statistics (mean, standard deviation) per month per site
long.data %>%
  group_by(Month, Site) %>%
  get_summary_stats(Stings, type = "mean_sd")


# order months in consecutive fashion, not alphabetically
data.sub$Month <- factor(data.sub$Month , levels=c("MAY_2021", "JULY_2021", "AUGUST_2021", "OCTOBER_2021",
                                                "NOVEMBER_2021"))


# plot boxplots over month, this plots combined sites per month 
# how to separate by site? 
boxplot.site <- ggplot(data.sub, aes(Month, Stings.on.Flag.cont, fill = Site)) +
  geom_boxplot() + scale_fill_manual(values=c("cornflowerblue", "deeppink"))

boxplot.site
# calculate means per month per site 
means <- aggregate(Stings.on.Flag.cont ~ Month*Site, data.sub, mean)


# plot means on top of boxplot 
boxplot.site.figure <- ggplot(data=data.sub, aes(x=Month, y=Stings.on.Flag.cont, fill=Site)) + geom_boxplot() +
  stat_summary(fun=mean, colour="darkred", geom="point", 
               shape=18, size=3, show.legend=FALSE, position = position_dodge2(width = 0.75, preserve = "single")) 
  geom_text(data = means, aes(label = Stings.on.Flag.cont, y = Stings.on.Flag.cont + 0.08)) 

# erase grey background and add title and custom color palette 
boxplot.site.figure <- boxplot.site.figure + theme_bw() +
  ggtitle("Stings on Flag Across Month and Between Sites with Plotted Average") +
  scale_fill_manual(values=c("cornflowerblue", "deeppink"))

boxplot.site.figure

## export figure
# Set the high-resolution png file 
png(("figures/Boxplot.SOF.Across.Month.and.Between.Sites.averages.jpg"),
    width = 10, 
    height = 7, 
    res = 300, 
    units = 'in')

#_______________________________________________________________________

# Biostatistical Analysis 

# Repeated Measures ANOVA to assess difference in mean between sites 
# across repeated measures through time. 

# Hypothesis: Honey bee colonies exhibit heightened defensiveness with \
# increasing time in ECR but not BFS. 

# assumptions of repeated measure ANOVA: normality, sphericity 



## Begin with a two-way repeated measures ANOVA to investigate whether there is \
# a significant interaction between site and month on DRV behavior score. 

# Compute the repeated measures ANOVA 
# dv = dependent variable 
# wid = unique individual identifier 
# within = within-subjects factor variables
# between = between-subjects factor variables
two.way.rep.measures <- anova_test(data = long.data,
                                   dv = Stings, 
                                   wid = Colony,
                                   within = Month, 
                                   between = Site)

get_anova_table(two.way.rep.measures)
# Site: p=0.000321 
# Month: p=0.007000
# Site:Month: p=0.01000


# Effect of treatment at each time point
one.way <- long.data %>%
  group_by(Month) %>%
  anova_test(dv = Stings, wid = Colony, between = Site) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# Site was significant in May, not in July or August, and then significant \
# in October and November 


# post hoc pairwise comparisons between treatment groups
pwc <- long.data %>%
  group_by(Month) %>%
  pairwise_t_test(
    Stings ~ Site, paired = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc

# BFS and ECR were significantly different in May, not in July or August, \
# and then significant in October and November 

# Effect of time at each level of treatment
one.way2 <- long.data %>%
  group_by(Site) %>%
  anova_test(dv = Stings, wid = Colony, within = Month) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

# Pairwise comparisons between time points
pwc2 <- long.data %>%
  group_by(Site) %>%
  pairwise_t_test(
    Stings ~ Month, paired = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc2



#_______________________________________________________________________



# More biostatistical analysis
# this time, we're running a generalized linear mixed model 

# the model 
# data in long format 
lme.sq <-lme(Stings ~ quality, random =~1|plot, data=long.data)




install.packages("lmerTest")

mix = lmer(score ~ Site + (1 | Colony) + time, data = long.data)

library()
rand(mix)
summary(mix)
