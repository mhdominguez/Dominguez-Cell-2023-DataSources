 
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
resultTable <- read.csv(file = 'JCF-FHF-SHF.tsv', sep='\t')
resultTable$Begin.Y.Rev = -resultTable$Begin.Y

#no re-scale needed, SVF exports in micron units
#resultTable$Begin.X = resultTable[,c("Begin.X")] / MaMuT_SVF_scale_factor
#resultTable$Begin.Y = resultTable[,c("Begin.Y")] / MaMuT_SVF_scale_factor
#resultTable$Begin.Z = resultTable[,c("Begin.Z")] / MaMuT_SVF_scale_factor
#resultTable$End.X = resultTable[,c("End.X")] / MaMuT_SVF_scale_factor
#resultTable$End.Y = resultTable[,c("End.Y")] / MaMuT_SVF_scale_factor
#resultTable$End.Z = resultTable[,c("End.Z")] / MaMuT_SVF_scale_factor

#==============
#look at unshifted data for anterior third of cells
#==============
colorpalette <- c("#F98A81","#00596f","#b5a96b")
resultTable$Region <- factor(resultTable$Region, levels=c("JCF", "aFHF", "aSHF"))

#Plot 2hr MA velocity
x_labels <- c("JCF", "aFHF", "aSHF") #unique(resultTable$Region)
resultTable %>%
filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
#filter( End.T > 60 ) %>%
filter( Begin.Y < 292 ) %>% #> 292 & Begin.Y < 442 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Region, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color=Region)) +
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
  scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  scale_y_continuous(limits=c(0,17.5)) +
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

test_1 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y < 292)
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y < 292 )
test_3  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y < 292 )
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_1$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)

#Plot displacement
x_labels <- c("JCF", "aFHF", "aSHF") #unique(resultTable$Region)
resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  filter( End.T > 60 ) %>%
  filter( Begin.Y < 292 ) %>% #> 292 & Begin.Y < 442 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Region, y=Peak.Displacement, alpha=0.6, color=Region)) +
  geom_violin(trim=FALSE,aes(fill = Region),color="black",size=1.5) +
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
  scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  scale_y_continuous(limits=c(0,150)) +
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

test_1 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y < 292)
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y < 292 )
test_3  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y < 292 )
t.test(test_1$Peak.Displacement, test_2$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Peak.Displacement, test_3$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_1$Peak.Displacement, test_3$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Peak.Displacement, test_2$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Peak.Displacement, test_3$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_1$Peak.Displacement, test_3$Peak.Displacement,  alternative = "two.sided", paired = FALSE)




#==============
#look at unshifted data for posterior two thirds of cells
#==============
colorpalette <- c("#F98A81","#00596f","#b5a96b")
resultTable$Region <- factor(resultTable$Region, levels=c("JCF", "aFHF", "aSHF"))

#Plot 2hr MA velocity
x_labels <- c("JCF", "aFHF", "aSHF") #unique(resultTable$Region)
resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  #filter( End.T > 60 ) %>%
  filter( Begin.Y > 292 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Region, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color=Region)) +
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
  scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3,scales="free_x") +
  #scale_y_continuous(trans='sqrt') +
  scale_y_continuous(limits=c(0,17.5)) +
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

test_1 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y > 292)
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y > 292 )
test_3  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y > 292 )
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_1$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)

#Plot displacement
x_labels <- c("JCF", "aFHF", "aSHF") #unique(resultTable$Region)
resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  filter( End.T > 60 ) %>%
  filter( Begin.Y > 292 ) %>% #> 292 & Begin.Y < 442 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Region, y=Peak.Displacement, alpha=0.6, color=Region)) +
  geom_violin(trim=FALSE,aes(fill = Region),color="black",size=1.5) +
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
  scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  scale_y_continuous(limits=c(0,150)) +
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

test_1 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y > 292)
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y > 292 )
test_3  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y > 292 )
t.test(test_1$Peak.Displacement, test_2$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Peak.Displacement, test_3$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_1$Peak.Displacement, test_3$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Peak.Displacement, test_2$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Peak.Displacement, test_3$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_1$Peak.Displacement, test_3$Peak.Displacement,  alternative = "two.sided", paired = FALSE)


#==============
###testing anterior to posterior
#==============
test_1 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y < 292)
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y < 292 )
test_3  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y < 292 )
test_4 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y > 292)
test_5  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y > 292 )
test_6  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y > 292 )
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_4$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Avg.Sliding.Velocity..micron.hr., test_5$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_3$Avg.Sliding.Velocity..micron.hr., test_6$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Avg.Sliding.Velocity..micron.hr., test_4$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Avg.Sliding.Velocity..micron.hr., test_5$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_3$Avg.Sliding.Velocity..micron.hr., test_6$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)

t.test(test_1$Peak.Displacement, test_4$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Peak.Displacement, test_5$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_3$Peak.Displacement, test_6$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Peak.Displacement, test_4$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Peak.Displacement, test_5$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_3$Peak.Displacement, test_6$Peak.Displacement,  alternative = "two.sided", paired = FALSE)



#==============
#Plot displacement and velocity vs. y-pos on scatterplot
#==============
x_labels <- c("JCF", "aFHF", "aSHF") #unique(resultTable$Region)

resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  #filter( End.T > 60 ) %>%
  #filter( Begin.Y < 292 ) %>% #> 292 & Begin.Y < 442 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Begin.Y.Rev, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color=Region)) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_fill_manual( values=colorpalette) +
  scale_color_manual(values=colorpalette) +
  geom_jitter(aes(color=Region), size=4, alpha=0.6) +
  geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=1.25) +
  #geom_crossbar(color="black", stat = "summary", fun=mean,size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(n.breaks=8) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
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
    panel.spacing = unit(60, "points"),
    panel.grid.major.x = element_line(colour="#DDDDDD",size=3),
    panel.grid.minor.x = element_line(colour="#DDDDDD",size=3),
    rect = element_blank(),
    plot.margin = margin(0, 30, 0, 0)
    #rect = element_rect( color="white", size=1)   
  ) +  
  xlab("")



resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  #filter( End.T > 60 ) %>%
  #filter( Begin.Y < 292 ) %>% #> 292 & Begin.Y < 442 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Begin.Y.Rev, y=Peak.Displacement, alpha=0.6, color=Region)) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_fill_manual( values=colorpalette) +
  scale_color_manual(values=colorpalette) +
  geom_jitter(aes(color=Region), size=4, alpha=0.6) +
  geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=1.25) +
  #geom_crossbar(color="black", stat = "summary", fun=mean,size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(n.breaks=8) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
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
    panel.spacing = unit(60, "points"),
    panel.grid.major.x = element_line(colour="#DDDDDD",size=3),
    panel.grid.minor.x = element_line(colour="#DDDDDD",size=3),
    rect = element_blank(),
    plot.margin = margin(0, 30, 0, 0)
    #rect = element_rect( color="white", size=1)   
  ) +  
  xlab("")



#==============
####Now look at shifted data
#==============
resultTable <- read.csv(file = 'JCF-FHF-SHF corrected.tsv', sep='\t')
resultTable$Begin.Y.Rev = -resultTable$Begin.Y
colorpalette <- c("#F98A81","#00596f","#b5a96b")
resultTable$Region <- factor(resultTable$Region, levels=c("JCF", "aFHF", "aSHF"))

#==============
#look at shifted data for anterior third of cells
#==============


#Plot 2hr MA velocity
x_labels <- c("JCF", "aFHF", "aSHF") #unique(resultTable$Region)
resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  #filter( End.T > 60 ) %>%
  filter( Begin.Y < 292 ) %>% #> 292 & Begin.Y < 442 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Region, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color=Region)) +
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
  scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  scale_y_continuous(limits=c(0,17.5)) +
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

test_1 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y < 292)
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y < 292 )
test_3  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y < 292 )
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_1$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)

#Plot displacement
x_labels <- c("JCF", "aFHF", "aSHF") #unique(resultTable$Region)
resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  filter( End.T > 60 ) %>%
  filter( Begin.Y < 292 ) %>% #> 292 & Begin.Y < 442 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Region, y=Peak.Displacement, alpha=0.6, color=Region)) +
  geom_violin(trim=FALSE,aes(fill = Region),color="black",size=1.5) +
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
  scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  scale_y_continuous(limits=c(0,150)) +
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

test_1 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y < 292)
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y < 292 )
test_3  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y < 292 )
t.test(test_1$Peak.Displacement, test_2$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Peak.Displacement, test_3$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_1$Peak.Displacement, test_3$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Peak.Displacement, test_2$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Peak.Displacement, test_3$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_1$Peak.Displacement, test_3$Peak.Displacement,  alternative = "two.sided", paired = FALSE)




#==============
#look at shifted data for posterior two thirds of cells
#==============
colorpalette <- c("#F98A81","#00596f","#b5a96b")
resultTable$Region <- factor(resultTable$Region, levels=c("JCF", "aFHF", "aSHF"))

#Plot 2hr MA velocity
x_labels <- c("JCF", "aFHF", "aSHF") #unique(resultTable$Region)
resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  #filter( End.T > 60 ) %>%
  filter( Begin.Y > 292 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Region, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color=Region)) +
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
  scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3,scales="free_x") +
  #scale_y_continuous(trans='sqrt') +
  scale_y_continuous(limits=c(0,17.5)) +
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

test_1 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y > 292)
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y > 292 )
test_3  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y > 292 )
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Avg.Sliding.Velocity..micron.hr., test_2$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_1$Avg.Sliding.Velocity..micron.hr., test_3$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)

#Plot displacement
x_labels <- c("JCF", "aFHF", "aSHF") #unique(resultTable$Region)
resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  filter( End.T > 60 ) %>%
  filter( Begin.Y > 292 ) %>% #> 292 & Begin.Y < 442 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Region, y=Peak.Displacement, alpha=0.6, color=Region)) +
  geom_violin(trim=FALSE,aes(fill = Region),color="black",size=1.5) +
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
  scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  scale_y_continuous(limits=c(0,150)) +
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

test_1 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y > 292)
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y > 292 )
test_3  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y > 292 )
t.test(test_1$Peak.Displacement, test_2$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Peak.Displacement, test_3$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_1$Peak.Displacement, test_3$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Peak.Displacement, test_2$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Peak.Displacement, test_3$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_1$Peak.Displacement, test_3$Peak.Displacement,  alternative = "two.sided", paired = FALSE)



###testing anterior to posterior
test_1 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y < 292)
test_2  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y < 292 )
test_3  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y < 292 )
test_4 <- filter(resultTable, str_detect(Region,"\\bJCF\\b") & End.T > 60 & Begin.Y > 292)
test_5  <- filter(resultTable, str_detect(Region,"aFHF") & End.T > 60 & Begin.Y > 292 )
test_6  <- filter(resultTable, str_detect(Region,"aSHF") & End.T > 60 & Begin.Y > 292 )
t.test(test_1$Avg.Sliding.Velocity..micron.hr., test_4$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Avg.Sliding.Velocity..micron.hr., test_5$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_3$Avg.Sliding.Velocity..micron.hr., test_6$Avg.Sliding.Velocity..micron.hr., alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Avg.Sliding.Velocity..micron.hr., test_4$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Avg.Sliding.Velocity..micron.hr., test_5$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)
var.test(test_3$Avg.Sliding.Velocity..micron.hr., test_6$Avg.Sliding.Velocity..micron.hr.,  alternative = "two.sided", paired = FALSE)

t.test(test_1$Peak.Displacement, test_4$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_2$Peak.Displacement, test_5$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
t.test(test_3$Peak.Displacement, test_6$Peak.Displacement, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
var.test(test_1$Peak.Displacement, test_4$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_2$Peak.Displacement, test_5$Peak.Displacement,  alternative = "two.sided", paired = FALSE)
var.test(test_3$Peak.Displacement, test_6$Peak.Displacement,  alternative = "two.sided", paired = FALSE)


#==============
#Plot displacement and velocity vs. y-pos on scatterplot
#==============
x_labels <- c("JCF", "aFHF", "aSHF") #unique(resultTable$Region)
resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  #filter( End.T > 60 ) %>%
  #filter( Begin.Y < 292 ) %>% #> 292 & Begin.Y < 442 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Begin.Y.Rev, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color=Region)) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_fill_manual( values=colorpalette) +
  scale_color_manual(values=colorpalette) +
  geom_jitter(aes(color=Region), size=4, alpha=0.6) +
  geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=1.25) +
  #geom_crossbar(color="black", stat = "summary", fun=mean,size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(n.breaks=8) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
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
    panel.spacing = unit(60, "points"),
    panel.grid.major.x = element_line(colour="#DDDDDD",size=3),
    panel.grid.minor.x = element_line(colour="#DDDDDD",size=3),
    rect = element_blank(),
    plot.margin = margin(0, 30, 0, 0)
    #rect = element_rect( color="white", size=1)   
  ) +  
  xlab("")


resultTable %>%
  filter( str_detect(Region,"\\bJCF\\b") | str_detect(Region,"aFHF") | str_detect(Region,"aSHF") ) %>%
  #filter( End.T > 60 ) %>%
  #filter( Begin.Y < 292 ) %>% #> 292 & Begin.Y < 442 ) %>% #filter out top third cells of crescent
  ggplot( aes(x=Begin.Y.Rev, y=Peak.Displacement, alpha=0.6, color=Region)) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_fill_manual( values=colorpalette) +
  scale_color_manual(values=colorpalette) +
  geom_jitter(aes(color=Region), size=4, alpha=0.6) +
  geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=1.25) +
  #geom_crossbar(color="black", stat = "summary", fun=mean,size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(n.breaks=8) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
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
    panel.spacing = unit(60, "points"),
    panel.grid.major.x = element_line(colour="#DDDDDD",size=3),
    panel.grid.minor.x = element_line(colour="#DDDDDD",size=3),
    rect = element_blank(),
    plot.margin = margin(0, 30, 0, 0)
    #rect = element_rect( color="white", size=1)   
  ) +  
  xlab("")

