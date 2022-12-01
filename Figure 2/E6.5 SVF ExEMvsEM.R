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
#open E6.5 dataset for ExEM vs EM
#==============
resultTable <- read.csv(file = 'Result_ExEMvsEM.tsv', sep='\t')
resultTable$Region <- as.factor(resultTable$Region)
colorpalette <- c("#EE0000","#707070")

#==============
#generate bar graphs / violin plots
#==============
#Plot displacement by region
x_labels <- unique(resultTable$Region)
resultTable %>%
  #filter( End.X > Begin.X ) %>%
  ggplot( aes(x=Region, y=Avg.Sliding.Velocity..micron.hr., alpha=1, color=Region)) +
  geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
  scale_fill_manual(values=colorpalette) + 
  scale_color_manual(values=colorpalette) +   
  geom_jitter(aes(color=Region), size=4, alpha=0.9) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("",""))#,"","","","","","","")) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3,scales="free_x") +
  scale_y_continuous(trans='sqrt',limits=c(1,20)) +
  #scale_y_continuous(limits=c(2.5,12.5),trans='sqrt') +
  #scale_y_continuous(limits=c(20,120)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.x = element_text(size=32,color="black",family="Liberation Sans"),
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

#t-tests
test_1 <- filter(resultTable, str_detect(Region,"ExEM") )
test_2  <- filter(resultTable, str_detect(Region,"EM") )
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)




#==============
#open E6.5 dataset for different timepoints
#==============
resultTable <- read.csv(file = 'Result_EMtimepointSeries.tsv', sep='\t')
resultTable$Timepoint <- as.factor(resultTable$Timepoint)
colorpalette <- c("#b1364d","#d3642d","#fd7ca7","#ee0000")

#==============
#generate bar graphs / violin plots
#==============
#Plot displacement by time envelope
x_labels <- unique(resultTable$Timepoint)
resultTable %>%
  #filter( End.X > Begin.X ) %>%
  ggplot( aes(x=Timepoint, y=Avg.Sliding.Velocity..micron.hr., alpha=1, color=Timepoint)) +
  geom_violin(trim=TRUE,aes(fill = Timepoint),color="black",size=1.5) +
  scale_fill_manual(values=colorpalette) + 
  scale_color_manual(values=colorpalette) +   
  geom_jitter(aes(color=Timepoint), size=4, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  coord_flip() + 
  scale_x_discrete(limits = rev(x_labels)) + #, labels=c("",""))#,"","","","","","","")) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3,scales="free_x") +
  scale_y_continuous(trans='sqrt',limits=c(1,20)) +
  #scale_y_continuous(limits=c(2.5,12.5),trans='sqrt') +
  #scale_y_continuous(limits=c(20,120)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.x = element_text(size=32,color="black",family="Liberation Sans"),
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


#t-tests
test_1 <- filter(resultTable, str_detect(Timepoint,"00") )
test_2  <- filter(resultTable, str_detect(Timepoint,"50") )
test_3  <- filter(resultTable, str_detect(Timepoint,"100") )
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_3$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)

