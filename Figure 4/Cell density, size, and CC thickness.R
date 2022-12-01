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
microns_per_pixel = 0.380490284561
#==============
#open density data E7.25
#==============
summaryTable <-read_ods(path="cell density data.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)
summaryTable$Stage <- as.factor(summaryTable$Timepoint)

#==============
#cell density E7.25
#==============
x_labels <- c( "+0h", "+3h", "+6h", "+9h", "+12h")

summaryTable %>%
  filter( str_detect(Stage,"0") | str_detect(Stage,"3") | str_detect(Stage,"6") | str_detect(Stage,"9") | str_detect(Stage,"12") ) %>%
  ggplot( aes(x=Stage, y=Density, fill=Stage, color=Stage, alpha=0.4)) +
  #geom_boxplot( outlier.shape = NA, lwd=2 ) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=c("#000000")), size=10, alpha=0.8) +
  coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  scale_y_continuous(limits=c(15,20)) + #,expand = c(0, 0)) +
  scale_x_discrete(labels= rev(x_labels), limits=rev) + #,expand=c(0,1)) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  scale_color_manual(values=c("#808080", "#808080", "#808080", "#808080", "#808080", "#808080")) +
  #scale_color_manual(values=water) +
  #scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    plot.margin = margin(0, 20, 0, 10),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White")
  ) +  
  xlab("")

#==============
#open cell volume data E7.25
#==============
summaryTable <-read_ods(path="cell size data.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)
summaryTable$Stage <- as.factor(summaryTable$Timepoint)

#==============
#cell volume
#==============
x_labels <- c( "+0h", "+3h", "+6h", "+9h", "+12h")

summaryTable %>%
  filter( str_detect(Stage,"0") | str_detect(Stage,"3") | str_detect(Stage,"6") | str_detect(Stage,"9") | str_detect(Stage,"12") ) %>%
  ggplot( aes(x=Stage, y=AvgCellVolume, fill=Stage, color=Stage, alpha=0.4)) +
  #geom_boxplot( outlier.shape = NA, lwd=2 ) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=c("#000000")), size=10, alpha=0.8) +
  coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(limits=c(15,20)) + #,expand = c(0, 0)) +
  scale_x_discrete(labels= rev(x_labels), limits=rev) + #,expand=c(0,1)) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  scale_color_manual(values=c("#808080", "#808080", "#808080", "#808080", "#808080", "#808080")) +
  #scale_color_manual(values=water) +
  #scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    plot.margin = margin(0, 20, 0, 10),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White")
  ) +  
  xlab("")


#==============
#open CC thickness data E7.25
#==============
summaryTable <-read_ods(path="F6 crescent dorsal-ventral thickness.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)
summaryTable$Stage <- as.factor(summaryTable$Timepoint)
summaryTable$Depth <- summaryTable$Depth * microns_per_pixel
#==============
#CC thickness
#==============
x_labels <- c( "+0h", "+3h", "+6h", "+9h", "+12h")

summaryTable %>%
  filter( str_detect(Stage,"0") | str_detect(Stage,"3") | str_detect(Stage,"6") | str_detect(Stage,"9") | str_detect(Stage,"12") ) %>%
  ggplot( aes(x=Stage, y=Depth, fill=Stage, color=Stage, alpha=0.4)) +
  #geom_boxplot( outlier.shape = NA, lwd=2 ) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=c("#000000")), size=10, alpha=0.8) +
  coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  scale_y_continuous() + #,expand = c(0, 0)) +
  scale_x_discrete(labels= rev(x_labels), limits=rev) + #,expand=c(0,1)) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  scale_color_manual(values=c("#808080", "#808080", "#808080", "#808080", "#808080", "#808080")) +
  #scale_color_manual(values=water) +
  #scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    plot.margin = margin(0, 20, 0, 10),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White")
  ) +  
  xlab("")


