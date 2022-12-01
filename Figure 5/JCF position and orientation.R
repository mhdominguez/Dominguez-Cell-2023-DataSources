#install.packages('readODS')
#font_import()
library(ggplot2)
library(readODS)
library(tidyverse)
library(viridis)
library(dplyr)
library(extrafont)
library(ggridges)
library(circular)
loadfonts()

library("RColorBrewer")
display.brewer.all()
brewer.pal(n = 9, name = "Paired")

MaMuT_scale_factor = 3.942282
microns_per_pixel = 0.380490284561
#==============
#open density data E7.25
#==============
summaryTable <-read_ods(path="combined data CC and JCF.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)
summaryTable$Region <- as.factor(summaryTable$Region)
#summaryTable$AngleRaw <- summaryTable$`Angle Raw` + 180
summaryTable$PosXTrans <- abs(summaryTable$PosX - 0.5)
#summaryTable$Angle3 <- ( summaryTable$Angle3 +355 ) %% 360
summaryTable$Angle4 <- 90 - abs(summaryTable$Angle2 - 90)


#==============
#radial histogram of angles -- at the extreme sides of the crescent
#==============
summaryTable %>%
  #filter( PosX < 0.6 & PosX > 0.4 ) %>%
  filter( PosX > 0.9 | PosX < 0.1 ) %>%
  ggplot( aes(x=PosX, fill=Region, color=Region, alpha=0.8)) +
  geom_histogram(aes(x=Angle3,fill=Region),alpha=0.5,binwidth=10,size=2 ) + #,position = "stack") + 

#scale_x_continuous(breaks=seq(15,360,30),limits=c(0,360)) + 
  scale_x_continuous(limits=c(0,360),breaks=seq(5,185,30)) +
  scale_y_continuous(trans='sqrt') +
  coord_polar(theta="x", start=3*pi/2-0.0872665, direction=-1) + 
#theme_bw() +
  #geom_boxplot( outlier.shape = NA, lwd=2 ) +
  #geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
 # geom_jitter(aes(color=c("#000000")), size=10, alpha=0.8) +
#  coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
 # scale_y_continuous(limits=c(0,180)) + #,expand = c(0, 0)) +
 # scale_x_discrete(labels= rev(x_labels), limits=rev) + #,expand=c(0,1)) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  scale_color_manual(values=c("#a0a0a0", "#7b4a51" )) + #, "#808080", "#808080", "#808080", "#808080")) +
  #scale_color_manual(values=water) +
  scale_fill_manual(values=c("#a0a0a0", "#F98A81" )) + 
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    #axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=0.3),
    panel.grid.major.y = element_line(colour="#A0A0A0",size=0.8),
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

test_1 <- filter(summaryTable, str_detect(Region,"JCF"), PosX > 0.9 | PosX < 0.1 )
test_2 <- filter(summaryTable, str_detect(Region,"CC"), PosX > 0.9 | PosX < 0.1 )
test_1_circular <- as.circular(test_1$Angle3,units="degrees")
test_2_circular <- as.circular(test_2$Angle3,units="degrees")

watson.two.test(test_1_circular,test_2_circular) #,alpha=0.05)



#==============
#radial histogram of angles -- at the top arch of the crescent
#==============
summaryTable %>%
  filter( PosX < 0.6 & PosX > 0.4 ) %>%
  #filter( PosX > 0.9 | PosX < 0.1 ) %>%
  ggplot( aes(x=PosX, fill=Region, color=Region, alpha=0.8)) +
  geom_histogram(aes(x=Angle3,fill=Region),alpha=0.5,binwidth=10,size=2 ) + #,position = "stack") + 
  
  #scale_x_continuous(breaks=seq(15,360,30),limits=c(0,360)) + 
  scale_x_continuous(limits=c(0,360),breaks=seq(5,185,30)) +
  scale_y_continuous(trans='sqrt') +
  coord_polar(theta="x", start=3*pi/2-0.0872665, direction=-1) + 
  #theme_bw() +
  #geom_boxplot( outlier.shape = NA, lwd=2 ) +
  #geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  # geom_jitter(aes(color=c("#000000")), size=10, alpha=0.8) +
  #  coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  # scale_y_continuous(limits=c(0,180)) + #,expand = c(0, 0)) +
  # scale_x_discrete(labels= rev(x_labels), limits=rev) + #,expand=c(0,1)) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  scale_color_manual(values=c("#a0a0a0", "#7b4a51" )) + #, "#808080", "#808080", "#808080", "#808080")) +
  #scale_color_manual(values=water) +
  scale_fill_manual(values=c("#a0a0a0", "#F98A81" )) + 
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    #axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=0.3),
    panel.grid.major.y = element_line(colour="#A0A0A0",size=0.8),
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

test_1 <- filter(summaryTable, str_detect(Region,"JCF"), PosX < 0.6 & PosX > 0.4  )
test_2 <- filter(summaryTable, str_detect(Region,"CC"),  PosX < 0.6 & PosX > 0.4  )
test_1_circular <- as.circular(test_1$Angle3,units="degrees")
test_2_circular <- as.circular(test_2$Angle3,units="degrees")

watson.two.test(test_1_circular,test_2_circular) #,alpha=0.05)




#==============
#scatter plot of angles versus X position
#==============
summaryTable %>%
  #filter( PosX < 0.6 & PosX > 0.4 ) %>%
  #filter( PosX > 0.9 | PosX < 0.1 ) %>%
  ggplot( aes(x=PosXTrans, y=Angle4, fill=Region, color=Region, alpha=0.8)) +
  geom_point(aes(color=Region), size=4, alpha=0.6) + #,position = "stack") + 
  
  #scale_x_continuous(breaks=seq(15,360,30),limits=c(0,360)) + 
  #scale_x_continuous(limits=c(0,360),breaks=seq(5,185,30)) +
  scale_y_continuous(limits=c(0,90),breaks=seq(0,90,30)) +
  #coord_polar(theta="x", start=3*pi/2-0.0872665, direction=-1) + 
  #theme_bw() +
  #geom_boxplot( outlier.shape = NA, lwd=2 ) +
  #geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  # geom_jitter(aes(color=c("#000000")), size=10, alpha=0.8) +
  #  coord_flip() + 
  facet_wrap( ncol = 2, vars(Region) )+
  stat_cor(aes(label = ..rr.label..), label.y.npc = 0.9, label.x.npc = 0.5, geom="label", vjust=0.5,hjust=0.5, color = "black", size=19) +
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  # scale_y_continuous(limits=c(0,180)) + #,expand = c(0, 0)) +
  # scale_x_discrete(labels= rev(x_labels), limits=rev) + #,expand=c(0,1)) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  scale_color_manual(values=c("#a0a0a0", "#F98A81" )) + #, "#808080", "#808080", "#808080", "#808080")) +
  #scale_color_manual(values=water) +
  scale_fill_manual(values=c("#a0a0a0", "#F98A81" )) + 
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    #axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=0.3),
    panel.grid.major.y = element_line(colour="#A0A0A0",size=0.8),
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
#scatter plot of angles versus Y position
#==============
summaryTable %>%
  #filter( PosX < 0.6 & PosX > 0.4 ) %>%
  #filter( PosX > 0.9 | PosX < 0.1 ) %>%
  ggplot( aes(x=PosY, y=Angle4, fill=Region, color=Region, alpha=0.8)) +
  geom_point(aes(color=Region), size=4, alpha=0.6) + #,position = "stack") + 
  
  #scale_x_continuous(breaks=seq(15,360,30),limits=c(0,360)) + 
  #scale_x_continuous(limits=c(0,360),breaks=seq(5,185,30)) +
  scale_y_continuous(limits=c(0,90),breaks=seq(0,90,30)) +
  #coord_polar(theta="x", start=3*pi/2-0.0872665, direction=-1) + 
  #theme_bw() +
  #geom_boxplot( outlier.shape = NA, lwd=2 ) +
  #geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  # geom_jitter(aes(color=c("#000000")), size=10, alpha=0.8) +
  #  coord_flip() + 
  facet_wrap( ncol = 2, vars(Region) )+
  stat_cor(aes(label = ..rr.label..), label.y.npc = 0.9, label.x.npc = 0.8, geom="label", vjust=0.5,hjust=0.5, color = "black", size=19) +
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  # scale_y_continuous(limits=c(0,180)) + #,expand = c(0, 0)) +
  # scale_x_discrete(labels= rev(x_labels), limits=rev) + #,expand=c(0,1)) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  scale_color_manual(values=c("#a0a0a0", "#F98A81" )) + #, "#808080", "#808080", "#808080", "#808080")) +
  #scale_color_manual(values=water) +
  scale_fill_manual(values=c("#a0a0a0", "#F98A81" )) + 
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    #axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=0.3),
    panel.grid.major.y = element_line(colour="#A0A0A0",size=0.8),
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

