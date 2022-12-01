 
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
MaMuT_SVF_scale_factor = 1 #all SVF exports are in micron
MaMuT_TGMMvanilla_scale_factor = 3.942282


#==============
#inspect SVF track data
#==============
resultTable <- read.csv(file = '3-layer Start Position.tsv', sep='\t')

#no re-scale needed, SVF exports in micron units
#resultTable$Begin.X = resultTable[,c("Begin.X")] / MaMuT_SVF_scale_factor
#resultTable$Begin.Y = resultTable[,c("Begin.Y")] / MaMuT_SVF_scale_factor
#resultTable$Begin.Z = resultTable[,c("Begin.Z")] / MaMuT_SVF_scale_factor
#resultTable$End.X = resultTable[,c("End.X")] / MaMuT_SVF_scale_factor
#resultTable$End.Y = resultTable[,c("End.Y")] / MaMuT_SVF_scale_factor
#resultTable$End.Z = resultTable[,c("End.Z")] / MaMuT_SVF_scale_factor

#==============
#generate violin plots of begin Zcoordinate
#==============

colorpalette <- c("#00596f","#2ac800","#9d0a0a")

#Plot origin by Z coordinate for top of crescent
x_labels <- unique(resultTable$Region)
resultTable %>%
filter( str_detect(Region,"JCF_PD") | str_detect(Region,"aFHF") | str_detect(Region,"Endo") ) %>%
#filter( Begin.X< 320, Begin.X>280 ) %>%
filter( Begin.Y < 280 ) %>% #filter out top of crescent:  where Y is within 180-278 units from top, or in other words the top 100um of the crescent
  ggplot( aes(x=Region, y=Begin.Z, alpha=0.6, color=Region)) +
  geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_fill_manual( values=colorpalette) +
  scale_color_manual(values=colorpalette) +
  geom_jitter(aes(color=Region), size=4, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels), labels=c("","","","","","","","","")) +
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

test_1 <- filter(resultTable, str_detect(Region,"JCF_PD") & Begin.Y < 280 )
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & Begin.Y < 280  )
test_3  <- filter(resultTable, str_detect(Region,"Endo") & Begin.Y < 280  )
t.test(test_1$Begin.Z, test_2$Begin.Z, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Begin.Z, test_3$Begin.Z, alternative = "two.sided", paired = FALSE, var.equal = FALSE)



#Plot origin by Z coordinate for base of crescent
resultTable %>%
  filter( str_detect(Region,"JCF_PD") | str_detect(Region,"aFHF") | str_detect(Region,"Endo") ) %>%
  #filter( Begin.X< 320, Begin.X>280 ) %>%
  filter( Begin.Y>350 & Begin.Y<450 ) %>% #filter out top of crescent:  where Y is within 180-278 units from top, or in other words the top 25um of the crescent
  ggplot( aes(x=Region, y=Begin.Z, alpha=0.6, color=Region)) +
  geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_fill_manual( values=colorpalette) +
  scale_color_manual(values=colorpalette) +
  geom_jitter(aes(color=Region), size=4, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels), labels=c("","","","","","","","","")) +
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

test_1 <- filter(resultTable, str_detect(Region,"JCF_PD") & Begin.Y>350 & Begin.Y<450 )
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & Begin.Y>350 & Begin.Y<450  )
test_3  <- filter(resultTable, str_detect(Region,"Endo") & Begin.Y>350 & Begin.Y<450  )
t.test(test_1$Begin.Z, test_2$Begin.Z, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Begin.Z, test_3$Begin.Z, alternative = "two.sided", paired = FALSE, var.equal = FALSE)




