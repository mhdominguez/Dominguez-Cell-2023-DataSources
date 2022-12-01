#install.packages('readODS')
#font_import()
library(ggplot2)
library(readODS)
library(tidyverse)
library(viridis)
library(dplyr)
library(extrafont)
library(ggridges)
loadfonts()

library("RColorBrewer")
display.brewer.all()
brewer.pal(n = 9, name = "Paired")

#==============
#inspect 3x3 data
#==============
resultTable <- read.csv(file = 'E6.75 3x3.tsv', sep='\t')

#==============
#generate bar graphs / violin plots
#==============

#Plot birthdate by region
x_labels <- unique(resultTable$Region)
resultTable %>%
  ggplot( aes(x=Region, y=Begin.T, alpha=1, color=Region)) +
  geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
  scale_fill_brewer(palette="Paired") + 
  scale_color_brewer(palette="Paired") +   
  geom_jitter(aes(color=Region), size=4, alpha=0.9) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  coord_flip() + 
  scale_x_discrete(limits = rev(x_labels), labels=c("","","","","","","","","")) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3,scales="free_x") +
  #scale_y_continuous(trans='sqrt') +
  #scale_y_continuous(limits=c(2.5,12.5),trans='sqrt') +
  #scale_y_continuous(limits=c(20,120)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    #strip.text.x = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(colour="#DDDDDD",size=3),
    panel.grid.minor.x = element_line(colour="#DDDDDD",size=3),
    rect = element_blank(),
    plot.margin = margin(0, 30, 0, 0)
    #rect = element_rect( color="white", size=1)   
  ) +  
  xlab("")

#Plot max displacement by region
resultTable %>%
  ggplot( aes(x=Region, y=Avg.Sliding.Velocity..micron.hr., alpha=1, color=Region)) +
  geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
  #stat_boxplot() +
  scale_fill_brewer(palette="Paired") + 
  scale_color_brewer(palette="Paired") +   
  geom_jitter(aes(color=Region), size=4, alpha=0.9) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  coord_flip() + 
  scale_x_discrete(limits = rev(x_labels), labels=c("","","","","","","","","")) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3,scales="free_x") +
  #scale_y_continuous(trans='sqrt') +
  #scale_y_continuous(limits=c(2.5,12.5),trans='sqrt') +
  scale_y_continuous(limits=c(10,200)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    #strip.text.x = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(colour="#DDDDDD",size=3),
    panel.grid.minor.x = element_line(colour="#DDDDDD",size=3),
    plot.margin = margin(0, 30, 0, 0),
    rect = element_blank()
    #rect = element_rect( color="white", size=1)   
  ) +  
  xlab("")



#Plot density by region
x_labels <- unique(resultTable$Region)
resultTable %>%
  ggplot( aes(x=Region, y=resultTable$Avg.Density..cells.per.12.radii., alpha=1, color=Region)) +
  geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
  scale_fill_brewer(palette="Paired") + 
  scale_color_brewer(palette="Paired") +   
  geom_jitter(aes(color=Region), size=4, alpha=0.9) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  coord_flip() + 
  scale_x_discrete(limits = rev(x_labels), labels=c("","","","","","","","","")) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3,scales="free_x") +
  #scale_y_continuous(trans='sqrt') +
  #scale_y_continuous(limits=c(2.5,12.5),trans='sqrt') +
  #scale_y_continuous(limits=c(20,120)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    #axis.text.x = element_text(size=36,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    #strip.text.x = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(colour="#DDDDDD",size=3),
    panel.grid.minor.x = element_line(colour="#DDDDDD",size=3),
    rect = element_blank()
    #rect = element_rect( color="white", size=1)   
  ) +  
  xlab("")




#t-test
Pos11 <- resultTable[which(resultTable$Region=="1-1"),]
Pos12 <- resultTable[which(resultTable$Region=="1-2"),]
Pos13 <- resultTable[which(resultTable$Region=="1-3"),]
Pos21 <- resultTable[which(resultTable$Region=="2-1"),]
Pos22 <- resultTable[which(resultTable$Region=="2-2"),]
Pos23 <- resultTable[which(resultTable$Region=="2-3"),]
Pos31 <- resultTable[which(resultTable$Region=="3-1"),]
Pos32 <- resultTable[which(resultTable$Region=="3-2"),]
Pos33 <- resultTable[which(resultTable$Region=="3-3"),]

Pos1X <- resultTable[grep("1-",resultTable$Region), ]
Pos2X <- resultTable[grep("2-",resultTable$Region), ]
Pos3X <- resultTable[grep("3-",resultTable$Region), ]
PosX1 <- resultTable[grep("-1",resultTable$Region), ]
PosX2 <- resultTable[grep("-2",resultTable$Region), ]
PosX3 <- resultTable[grep("-3",resultTable$Region), ]

t.test(Pos1X$Begin.T,Pos3X$Begin.T)
t.test(PosX1$Begin.T,PosX3$Begin.T)

#t.test(Pos21$Begin.T,Pos23$Begin.T)
#t.test(Pos31$Begin.T,Pos33$Begin.T)
#t.test(Pos21$Begin.T,Pos23$Peak.Displacement)
#t.test(Pos31$Begin.T,Pos33$Peak.Displacement)

t.test(Pos1X$Avg.Density..cells.per.12.radii.,Pos3X$Avg.Density..cells.per.12.radii.)
t.test(PosX1$Avg.Density..cells.per.12.radii.,PosX3$Avg.Density..cells.per.12.radii.)


t.test(Pos1X$Avg.Sliding.Velocity..micron.hr.,Pos3X$Avg.Sliding.Velocity..micron.hr.)
t.test(PosX1$Avg.Sliding.Velocity..micron.hr.,PosX3$Avg.Sliding.Velocity..micron.hr.)
#t.test(Pos11$Begin.T,Pos31$Begin.T)
#t.test(Pos13$Begin.T,Pos33$Begin.T)
#t.test(Pos12$Begin.T,Pos32$Begin.T)


