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

#open data
summaryTable <-read_ods(path="cell counts data.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)
summaryTable$Stage <- as.factor(summaryTable$Timepoint)

#==============
#cell count
#==============
x_labels <- c( "E6.5", "E6.75", "E7.0", "E7.25", "E7.5", "E7.75", "E8.0")
new_xy_labels <- c("x" = "Anterior-Posterior", "y" = "Medial-Lateral")
new_13_labels <- c("1" = "", "3" = "")
fills65 <-c("#64ABE3","#FDD8B5")
fills67 <-c("#15B2D1","#DFCE9D")
fills70 <-c("#65CBDA","#F9D199")
sand <- c(fills65[2],fills67[2],fills70[2])
water <- c(fills65[1],fills67[1],fills70[1])

summaryTable %>%
  ggplot( aes(x=Stage, y=Mesp1, fill=Stage, color=Stage, alpha=0.4)) +
  #geom_boxplot( outlier.shape = NA, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=c("#a00050")), size=10, alpha=0.8) +
  #coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  scale_y_continuous(trans="log") + #,expand = c(0, 0)) +
  scale_x_discrete(labels= x_labels) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  scale_color_manual(values=c("#a00050", "#a00050", "#a00050", "#a00050", "#a00050", "#a00050", "#a00050")) +
  #scale_color_manual(values=water) +
  #scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
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


summaryTable %>%
  ggplot( aes(x=Stage, y=F6, fill=Stage, color=Stage, alpha=0.4)) +
  #geom_boxplot( outlier.shape = NA, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=c("#0fa85c")), size=10, alpha=0.8) +
  #coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  scale_y_continuous(trans="log") + #,expand = c(0, 0)) +
  scale_x_discrete(labels= x_labels) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  scale_color_manual(values=c("#0fa85c", "#0fa85c", "#0fa85c", "#0fa85c", "#0fa85c", "#0fa85c", "#0fa85c")) +
  #scale_color_manual(values=water) +
  #scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
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


