#! /usr/local/bin/R

setwd("/Users/tpd0001/Desktop/R_tut")


liver <- read.csv("r_tut.liver_drug.txt")

# Student's t test for drug vs no_drug in liver
t.test(count ~ treatment, data=liver)



spleen <- read.csv("r_tut.spleen_drug.txt")

# one-way ANOVA
aov(count ~ treatment, data = spleen)
# or 
res.aov <- aov(count ~ treatment, data = spleen)
summary(res.aov)

# two-way ANOVA
aov(count ~ gene * treatment, data = spleen)
or 
res.aov2 <- aov(count ~ gene * treatment, data = spleen)
summary(res.aov2)

# Tukey multiple comparsions of means
TukeyHSD(res.aov2, which = "gene")

