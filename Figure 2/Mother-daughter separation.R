#install.packages('readODS')
#font_import()
library(ggplot2)
library(readODS)
library(tidyverse)
library(viridis)
library(dplyr)
library(extrafont)
library(ggridges)
library(ggpubr)
loadfonts()

library("RColorBrewer")
display.brewer.all()
brewer.pal(n = 9, name = "Paired")

MaMuT_scale_factor = 3.942282

#open E6.5 dataset
resultTable65 <- read.csv(file = 'E6.5_daughter_separation_data.tsv', sep='\t')
resultTable65$Stage = "E6.5"
resultTable65$IdxDDAP = resultTable65[,c("Daughter.Daugther.Distance")] * 100 / ( 226.83 * MaMuT_scale_factor )
resultTable65$IdxMDVel = resultTable65[,c("Mother.Daughter.Distance")] * 100 / ( 2 * 35.96 ) #no need to scale since everything in MaMuT terms
resultTable65$IdxDDVel = resultTable65[,c("Daughter.Daugther.Distance")] * 100 / ( 2 * 35.96 ) #no need to scale since everything in MaMuT terms


resultTable65$Daughter.Daugther.Distance = resultTable65[,c("Daughter.Daugther.Distance")] / MaMuT_scale_factor

#open E6.75 dataset
resultTable67 <- read.csv(file = 'E6.75_daughter_separation_data.tsv', sep='\t')
resultTable67$Stage = "E6.75"
resultTable67$IdxDDAP = resultTable67[,c("Daughter.Daugther.Distance")] * 100 / ( 265.5 * MaMuT_scale_factor )
resultTable67$IdxMDVel = resultTable67[,c("Mother.Daughter.Distance")] * 100 / ( 2 * 33.71 ) #no need to scale since everything in MaMuT terms
resultTable67$IdxDDVel = resultTable67[,c("Daughter.Daugther.Distance")] * 100 / ( 2 * 33.71 ) #no need to scale since everything in MaMuT terms


resultTable67$Daughter.Daugther.Distance = resultTable67[,c("Daughter.Daugther.Distance")] / MaMuT_scale_factor

#open E7.0 dataset
resultTable70 <- read.csv(file = 'E7.0_daughter_separation_data.tsv', sep='\t')
resultTable70$Stage = "E7.0"
resultTable70$IdxDDAP = resultTable70[,c("Daughter.Daugther.Distance")] * 100 / ( 400.3 * MaMuT_scale_factor )
resultTable70$IdxMDVel = resultTable70[,c("Mother.Daughter.Distance")] * 100 / ( 2 * 35.54 ) #no need to scale since everything in MaMuT terms
resultTable70$IdxDDVel = resultTable70[,c("Daughter.Daugther.Distance")] * 100 / ( 2 * 35.54 ) #no need to scale since everything in MaMuT terms


resultTable70$Daughter.Daugther.Distance = resultTable70[,c("Daughter.Daugther.Distance")] / MaMuT_scale_factor

#merge datasets
resultTable <- rbind(resultTable65,resultTable67,resultTable70 )
resultTable$Stage <- as.factor(resultTable$Stage)

#labels and colors
x_labels <- c( "E6.5", "E6.75", "E7.0")
fills65 <-c("#64ABE3","#FDD8B5")
fills67 <-c("#15B2D1","#DFCE9D")
fills70 <-c("#65CBDA","#F9D199")
sand <- c(fills65[2],fills67[2],fills70[2])
water <- c(fills65[1],fills67[1],fills70[1])

#not used here
new_xy_labels <- c("x" = "Anterior-Posterior", "y" = "Medial-Lateral")
new_13_labels <- c("1" = "", "3" = "")

#==============
#summary table plots
#==============

resultTable %>%
  ggplot( aes(x=Stage, y=Daughter.Daugther.Distance, fill=Stage, color=Stage, alpha=0.4)) +
  geom_violin( trim=FALSE, outlier.shape = NA, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels, limits = rev) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=water) +
  scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White")
  ) +  
  xlab("")

resultTable %>%
  ggplot( aes(x=Stage, y=Mother.Daughter.Distance, fill=Stage, color=Stage, alpha=0.4)) +
  geom_violin( trim=FALSE, outlier.shape = NA, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels, limits = rev) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=water) +
  scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White")
  ) +  
  xlab("")

resultTable %>%
  ggplot( aes(x=Stage, y=IdxMDVel, fill=Stage, color=Stage, alpha=0.4)) +
  geom_violin( trim=FALSE, outlier.shape = NA, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels, limits = rev) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=water) +
  scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White")
  ) +  
  xlab("")

resultTable %>%
  ggplot( aes(x=Stage, y=IdxDDAP, fill=Stage, color=Stage, alpha=0.4)) +
  geom_violin( trim=FALSE, outlier.shape = NA, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels, limits = rev) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=water) +
  scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White")
  ) +  
  xlab("")

resultTable %>%
  ggplot( aes(x=Stage, y=IdxDDVel, fill=Stage, color=Stage, alpha=0.4)) +
  geom_violin( trim=FALSE, outlier.shape = NA, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels, limits = rev) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=water) +
  scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White")
  ) +  
  xlab("")



#==============
#t-test
#==============
PosE65 <- resultTable[which(resultTable$Stage=="E6.5"),]
PosE67 <- resultTable[which(resultTable$Stage=="E6.75"),]
PosE70 <- resultTable[which(resultTable$Stage=="E7.0"),]

t.test(PosE65$Mother.Daughter.Distance,PosE67$Mother.Daughter.Distance)
t.test(PosE67$Mother.Daughter.Distance,PosE70$Mother.Daughter.Distance)

t.test(PosE65$Daughter.Daugther.Distance,PosE67$Daughter.Daugther.Distance)
t.test(PosE67$Daughter.Daugther.Distance,PosE70$Daughter.Daugther.Distance)

t.test(PosE65$IdxDDAP,PosE67$IdxDDAP)
t.test(PosE67$IdxDDAP,PosE70$IdxDDAP)

t.test(PosE65$IdxMDVel,PosE67$IdxMDVel)
t.test(PosE67$IdxMDVel,PosE70$IdxMDVel)

t.test(PosE65$IdxDDVel,PosE67$IdxDDVel)
t.test(PosE67$IdxDDVel,PosE70$IdxDDVel)
