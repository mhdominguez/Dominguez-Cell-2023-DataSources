#install.packages('readODS')
#font_import()
library(ggplot2)
library(readODS)
library(tidyverse)
library(viridis)
library(dplyr)
#library(extrafont)
library(patchwork)
library(reshape2)
#loadfonts()

library("RColorBrewer")
display.brewer.all()
brewer.pal(n = 9, name = "Paired")


#==============
#labels and colors
#==============
#paired_colors <- rev(brewer.pal(name = "Set1", n=6)[c(2,6)])
paired_colors <- c("Dark Gray","Purple")



#==============
# read data and add factors
#==============
resultTable_quant <-read_ods(path="FHF SHF control vs Isl1KO.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)
resultTable_quant$Condition <- as.factor(resultTable_quant$Condition)
resultTable_quant <- melt(resultTable_quant)

resultTable_count <- read.csv("Embryo Mef2c cell counts.csv")
resultTable_count$Condition <- as.factor(resultTable_count$Condition)






#==============
# pHH3 cell numbers
#==============
#Region.labels <- c("pHH3/FHF","pHH3/SHF")

p1 <- resultTable_quant %>% filter( variable == "pHH3 / SHF" | variable == "pHH3 / FHF" ) %>%
  ggplot( aes(x=Condition, y=value*100, color=Condition)) +
  #geom_violin( trim=FALSE, lwd=2 ) +
  #geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Condition), size=8, alpha=0.8, width=0.2) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25, show.legend=FALSE) +
  #coord_flip() + 
  facet_grid( ~variable,switch="both",labeller = labeller(Region = Region.labels))+
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete() + #labels= x_labels ) +#
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual (values=paired_colors) +
  #scale_fill_manual(values=paired_colors) +
  #guides(color = guide_legend(reverse=TRUE)) +
  theme(
    legend.position = "none",
    #legend.position = c(0.78, 0.21),
    legend.box.background = element_rect(colour="#A0A0A0",size=1),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
    legend.title = element_text(size=56,vjust=1,hjust=0.5,color="black"),
    legend.text = element_text(size=48,color="black", angle=0, vjust=0.65, hjust=0),
    plot.title = element_text(size=56,vjust=1,hjust=0.5,color="black", face="bold"),
    #axis.text.y = element_text(size=32,color="black"),
    axis.title = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black", face="bold"),
    axis.text.y = element_text(size=48,color="black", angle=0, vjust=0.65, hjust=0.5),
    axis.text.x = element_blank(), #element_text(size=60,color="black", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black", face="bold" ),
    strip.text = element_text(size=48,color="black"),
    #panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(6,"lines"),
    #panel.background = element_rect(fill="Light Gray"),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White"),
    strip.placement = "outside"
  ) +  
  xlab("")
ggsave( plot = p1, height=7, width=8.8, filename = "pHH3 quantifications.png", units = "in",dpi = 100 )

Isl1KO <- resultTable_quant %>% filter( Condition == "Isl1 KO" )
Control <- resultTable_quant %>% filter(  Condition == "Control" )

t.test(Isl1KO[which(Isl1KO$variable=="pHH3 / FHF"),]$value,   Control[which(Control$variable=="pHH3 / FHF"),]$value)
t.test(Isl1KO[which(Isl1KO$variable=="pHH3 / SHF"),]$value,   Control[which(Control$variable=="pHH3 / SHF"),]$value)




#==============
# Mef2c cell numbers
#==============
p2 <- resultTable_count %>% 
  ggplot( aes(x=Condition, y=Estimated.Mef2c.Count, color=Condition)) +
  #geom_violin( trim=FALSE, lwd=2 ) +
  #geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Condition), size=8, alpha=0.8, width=0.2) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25, show.legend=FALSE) +
  #coord_flip() + 
  #facet_grid( ~variable,switch="both",labeller = labeller(Region = Region.labels))+
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(0,5000)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete() + #labels= x_labels ) +#
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual (values=paired_colors) +
  #scale_fill_manual(values=paired_colors) +
  #guides(color = guide_legend(reverse=TRUE)) +
  theme(
    legend.position = "none",
    #legend.position = c(0.78, 0.21),
    legend.box.background = element_rect(colour="#A0A0A0",size=1),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
    legend.title = element_text(size=56,vjust=1,hjust=0.5,color="black"),
    legend.text = element_text(size=48,color="black", angle=0, vjust=0.65, hjust=0),
    plot.title = element_text(size=56,vjust=1,hjust=0.5,color="black", face="bold"),
    #axis.text.y = element_text(size=32,color="black"),
    axis.title = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black", face="bold"),
    axis.text.y = element_text(size=48,color="black", angle=0, vjust=0.65, hjust=0.5),
    axis.text.x = element_blank(), #element_text(size=60,color="black", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black", face="bold" ),
    strip.text = element_text(size=48,color="black"),
    #panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(6,"lines"),
    #panel.background = element_rect(fill="Light Gray"),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White"),
    strip.placement = "outside"
  ) +  
  xlab("")
ggsave( plot = p2, height=7, width=4.7, filename = "Mef2c quantifications.png", units = "in",dpi = 100 )

Isl1KO <- resultTable_count %>% filter( Condition == "Isl1 KO" )
Control <- resultTable_count %>% filter(  Condition == "Control" )

t.test(Isl1KO$Estimated.Mef2c.Count,   Control$Estimated.Mef2c.Count)



#==============
# FHF/SHF cell numbers
#==============
p3 <- resultTable_quant %>% filter( variable == "SHF / total" ) %>%
  ggplot( aes(x=Condition, y=value*100, color=Condition)) +
  #geom_violin( trim=FALSE, lwd=2 ) +
  #geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Condition), size=8, alpha=0.8, width=0.2) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25, show.legend=FALSE) +
  #coord_flip() + 
  #facet_grid( ~variable,switch="both",labeller = labeller(Region = Region.labels))+
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(30,50)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete() + #labels= x_labels ) +#
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual (values=paired_colors) +
  #scale_fill_manual(values=paired_colors) +
  #guides(color = guide_legend(reverse=TRUE)) +
  theme(
    legend.position = "none",
    #legend.position = c(0.78, 0.21),
    legend.box.background = element_rect(colour="#A0A0A0",size=1),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
    legend.title = element_text(size=56,vjust=1,hjust=0.5,color="black"),
    legend.text = element_text(size=48,color="black", angle=0, vjust=0.65, hjust=0),
    plot.title = element_text(size=56,vjust=1,hjust=0.5,color="black", face="bold"),
    #axis.text.y = element_text(size=32,color="black"),
    axis.title = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black", face="bold"),
    axis.text.y = element_text(size=48,color="black", angle=0, vjust=0.65, hjust=0.5),
    axis.text.x = element_blank(), #element_text(size=60,color="black", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black", face="bold" ),
    strip.text = element_text(size=48,color="black"),
    #panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(6,"lines"),
    #panel.background = element_rect(fill="Light Gray"),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White"),
    strip.placement = "outside"
  ) +  
  xlab("")
ggsave( plot = p3, height=7, width=4.7, filename = "SHF-ratio quantifications.png", units = "in",dpi = 100 )

Isl1KO <- resultTable_quant %>% filter( Condition == "Isl1 KO" )
Control <- resultTable_quant %>% filter(  Condition == "Control" )

t.test(Isl1KO[which(Isl1KO$variable=="SHF / total"),]$value,   Control[which(Control$variable=="SHF / total"),]$value)













