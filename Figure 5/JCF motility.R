#install.packages('readODS')
#font_import()
library(ggplot2)
library(readODS)
library(tidyverse)
library(viridis)
library(dplyr)
library(extrafont)
library(ggridges)
library(ggpbr)
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
summaryTable <-read_ods(path="JCF-CC 30min motility data.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)
summaryTable$Region <- as.factor(summaryTable$Region)
#summaryTable$AngleRaw <- summaryTable$`Angle Raw` + 180
#summaryTable$PosXTrans <- abs(summaryTable$PosX - 0.5)
#summaryTable$Angle3 <- ( summaryTable$Angle3 +355 ) %% 360
summaryTable$Avg.Sliding.Velocity..micron.hr. <- summaryTable$Displacement *2
summaryTable$Angle3 <- summaryTable$Theta + 180


#==============
#violin plot of motility
#==============

x_labels <- c("JCF", "CC") #unique(resultTable$Region)
summaryTable %>%
  #filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  #filter( Begin.X> 320, Begin.X<280 ) %>%
  #filter( End.T > 50 ) %>% #filter out top of crescent:  where Y is within 180-278 units from top, or in other words the top 100um of the crescent
  ggplot( aes(x=Region, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color=Region)) +
  geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_fill_manual( values=c("#a0a0a0", "#F98A81" )) +
  scale_color_manual(values=c("#a0a0a0", "#F98A81" )) +
  geom_jitter(aes(color=Region), size=4, alpha=0.6) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  coord_flip() + 
  scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
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

#t-tests
test_1 <- filter(summaryTable, str_detect(Region,"CC"))
test_2  <- filter(summaryTable,  str_detect(Region,"JCF") )
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)



#==============
#radial histogram of direction moved
#==============
summaryTable %>%
  #filter( XM < 300 & XM > 200 ) %>%
  #filter( PosX > 0.9 | PosX < 0.1 ) %>%
  ggplot( aes(x=Displacement, fill=Region, color=Region, alpha=0.8)) +
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

summaryTable %>% 
  group_by(Region) %>%
  summarise(no_rows = length(Region))


test_1 <- filter(summaryTable, str_detect(Region,"JCF"))
test_2 <- filter(summaryTable, str_detect(Region,"CC") )
test_1_circular <- as.circular(test_1$Angle3,units="degrees")
test_2_circular <- as.circular(test_2$Angle3,units="degrees")

watson.two.test(test_1_circular,test_2_circular) #,alpha=0.05)



