#! /usr/local/bin/R

library(ggbiplot)
library(ggplot2)
library(readxl)

setwd("~/Desktop/R_crash_course")


# Simple version
abdata <- read_excel("pca/pca_demo_min.xlsx")
abdata.pca <- prcomp(abdata, scale = TRUE, center = TRUE)
simple <- ggbiplot(abdata.pca)



# More complex version
abdata2 <- read_excel("pca/pca_demo.xlsx")
#abdata2.trees <- as.numeric(as.factor(as.character(abdata2$tree)))
abdata2.trees <- as.factor(as.character(abdata2$tree))
pcdata <- subset(abdata2, select = -c(1))

pcdata.pca <- prcomp(pcdata, scale = TRUE, center = TRUE)
grouped <- ggbiplot(pcdata.pca, groups = abdata2.trees)


mycolors <- c("#FF7E79", "#A5712D", "#7A81FF")
my_theme <- theme(legend.position = "none", panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
my_theme2 <- list(my_theme, scale_color_manual(values = mycolors))

glorious <- ggbiplot(pcdata.pca, obs.scale = 1, var.scale = 1, groups = abdata2.trees, ellipse = TRUE, var.axes = FALSE) + scale_color_manual(values = mycolors) + geom_point(aes(color = abdata2.trees), size = 3) + my_theme
