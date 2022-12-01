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

data_input <-read_ods(path="F6 lineage tracing quants.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)
data_input$Percent <- data_input$Percent * 100
data_input$Region <- as.factor(data_input$Region)
data_input$Stage <- as.factor(data_input$Stage)

#labels and colors
x_labels <- c( "E5.5", "E6.5", "E7.5")
fills65 <-c("#64ABE3","#FDD8B5")
fills67 <-c("#15B2D1","#DFCE9D")
fills70 <-c("#65CBDA","#F9D199")
sand <- c(fills65[2],fills67[2],fills70[2])
water <- c(fills65[1],fills67[1],fills70[1])
oocalc <- c("#004586","#FF420E","#FFD320","#579D1C","#7E0021","#83CAFF","#314004","#AECF00")

#==============
#Violin plots of all the stuff: raw MD (um separation after 2hr), raw DD (um separation after 2hr), indexd DD (% coverage of AP distance of embryo after 2hr), indexd MD (separation speed as % of average track speed), indexed DD2 (separation speed as % of average track speed)
#==============

#==============
#summary table plots
#==============
.f = function() {

data_input %>%
  ggplot( aes(x=Stage, y=Percent, alpha=0.4)) +
  #geom_boxplot( ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Region), size=12, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,aes(color=Region),size=2) +
  #coord_flip() + 
  facet_wrap( ncol = 4, vars(Region) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(0,42)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=oocalc) +
  scale_fill_manual(values=oocalc) +
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
    axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    panel.grid.major.x = element_blank(),
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

data_input %>%
filter( Region=="Atria" ) %>%
  ggplot( aes(x=Stage, y=Percent, alpha=0.4)) +
  #geom_boxplot( ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="#7E0021", size=14, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="#7E0021",size=3) +
  #coord_flip() + 
  #facet_wrap( ncol = 4, vars(Region) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(0,21)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  #scale_color_manual(values=oocalc) +
  #scale_fill_manual(values=oocalc) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=45, vjust=0.65),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    panel.grid.major.x = element_blank(),
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

data_input %>%
  filter( Region=="LV" ) %>%
  ggplot( aes(x=Stage, y=Percent, alpha=0.4)) +
  #geom_boxplot( ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="#004586", size=14, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="#004586",size=3) +
  #coord_flip() + 
  #facet_wrap( ncol = 4, vars(Region) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(0,42)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  #scale_color_manual(values=oocalc) +
  #scale_fill_manual(values=oocalc) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=45, vjust=0.65),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    panel.grid.major.x = element_blank(),
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

data_input %>%
  filter( Region=="RV" ) %>%
  ggplot( aes(x=Stage, y=Percent, alpha=0.4)) +
  #geom_boxplot( ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="#FF420E", size=14, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="#FF420E",size=3) +
  #coord_flip() + 
  #facet_wrap( ncol = 4, vars(Region) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(0,18)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  #scale_color_manual(values=oocalc) +
  #scale_fill_manual(values=oocalc) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=45, vjust=0.65),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    panel.grid.major.x = element_blank(),
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

data_input %>%
  filter( Region=="OFT" ) %>%
  ggplot( aes(x=Stage, y=Percent, alpha=0.4)) +
  #geom_boxplot( ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="#FFD320", size=14, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="#FFD320",size=3) +
  #coord_flip() + 
  #facet_wrap( ncol = 4, vars(Region) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(0,21)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  #scale_color_manual(values=oocalc) +
  #scale_fill_manual(values=oocalc) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=45, vjust=0.65),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=45, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    panel.grid.major.x = element_blank(),
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
}


data_input %>%
  filter( Region=="Atria" ) %>%
  ggplot( aes(x=Stage, y=Percent, alpha=0.4)) +
  #geom_boxplot( ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="#7E0021", size=14, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="#7E0021",size=3) +
  coord_flip() + 
  scale_x_discrete(limits = rev(x_labels)) +
  #facet_wrap( ncol = 4, vars(Region) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(0,21)) +
  #scale_y_continuous(trans="sqrt") +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  #scale_color_manual(values=oocalc) +
  #scale_fill_manual(values=oocalc) +
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
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    panel.grid.major.x = element_blank(),
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

data_input %>%
  filter( Region=="LV" ) %>%
  ggplot( aes(x=Stage, y=Percent, alpha=0.4)) +
  #geom_boxplot( ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="#004586", size=14, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="#004586",size=3) +
  coord_flip() + 
  scale_x_discrete(limits = rev(x_labels)) +
  #facet_wrap( ncol = 4, vars(Region) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(0,42)) +
  #scale_y_continuous(trans="sqrt") +
  #scale_x_discrete(labels= x_labels) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  #scale_color_manual(values=oocalc) +
  #scale_fill_manual(values=oocalc) +
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
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    panel.grid.major.x = element_blank(),
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

data_input %>%
  filter( Region=="RV" ) %>%
  ggplot( aes(x=Stage, y=Percent, alpha=0.4)) +
  #geom_boxplot( ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="#FF420E", size=14, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="#FF420E",size=3) +
  coord_flip() + 
  scale_x_discrete(limits = rev(x_labels)) +
  #facet_wrap( ncol = 4, vars(Region) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(0,18)) +
  #scale_y_continuous(trans="sqrt") +
  #scale_x_discrete(labels= x_labels) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  #scale_color_manual(values=oocalc) +
  #scale_fill_manual(values=oocalc) +
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
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    panel.grid.major.x = element_blank(),
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

data_input %>%
  filter( Region=="OFT" ) %>%
  ggplot( aes(x=Stage, y=Percent, alpha=0.4)) +
  #geom_boxplot( ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="#FFD320", size=14, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="#FFD320",size=3) +
  coord_flip() + 
  scale_x_discrete(limits = rev(x_labels)) +
  #facet_wrap( ncol = 4, vars(Region) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(0,21)) +
  #scale_y_continuous(trans="sqrt") +
  #scale_x_discrete(labels= x_labels) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  #scale_color_manual(values=oocalc) +
  #scale_fill_manual(values=oocalc) +
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
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=45, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    panel.grid.major.x = element_blank(),
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
#PosE65 <- resultTable[which(resultTable$Stage=="E6.5"),]
resultTable65 <- data_input[which(data_input$Stage=="E6.5"),]
resultTable75 <- data_input[which(data_input$Stage=="E7.5"),]

RV_E65 <- resultTable65[which(resultTable65$Region=="RV"),]
RV_E75 <- resultTable75[which(resultTable75$Region=="RV"),]

OFT_E65 <- resultTable65[which(resultTable65$Region=="OFT"),]
OFT_E75 <- resultTable75[which(resultTable75$Region=="OFT"),]

LV_E65 <- resultTable65[which(resultTable65$Region=="LV"),]
LV_E75 <- resultTable75[which(resultTable75$Region=="LV"),]

Atria_E65 <- resultTable65[which(resultTable65$Region=="Atria"),]
Atria_E75 <- resultTable75[which(resultTable75$Region=="Atria"),]

AVC_E65 <- resultTable65[which(resultTable65$Region=="AVC"),]
AVC_E75 <- resultTable75[which(resultTable75$Region=="AVC"),]

t.test(RV_E65$Percent,RV_E75$Percent)
t.test(LV_E65$Percent,LV_E75$Percent)
t.test(OFT_E65$Percent,OFT_E75$Percent)
t.test(Atria_E65$Percent,Atria_E75$Percent)
t.test(AVC_E65$Percent,AVC_E75$Percent)


