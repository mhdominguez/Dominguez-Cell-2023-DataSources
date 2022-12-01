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

untar("Track_crossing_x.tar.bz2")
untar("Track_crossing_y.tar.bz2")

#open E6.5 dataset
resultTable65x <- read.csv(file = 'E6.5 crossing x to t4.5 maxdist250.tsv', sep='\t')
resultTable65x[,'Region']<-factor(resultTable65x[,'Region'])
resultTable65x$Axis <- "x"
resultTable65y <- read.csv(file = 'E6.5 crossing y to t4.5 maxdist250.tsv', sep='\t')
resultTable65y[,'Region']<-factor(resultTable65y[,'Region'])
resultTable65y$Axis <- "y"
#resultTable65 <- dplyr::bind_rows(resultTable65x,resultTable65y)
resultTable65 <- rbind(resultTable65x,resultTable65y)
resultTable65$Begin.Diff = resultTable65[,c("Begin.Diff")] / MaMuT_scale_factor
resultTable65$End.Diff = resultTable65[,c("End.Diff")] / MaMuT_scale_factor
resultTable65$BE.Diff = resultTable65[,c("End.Diff")] - resultTable65[,c("Begin.Diff")]

#open E6.75 dataset
resultTable67x <- read.csv(file = 'E6.75 crossing x to t4.5 maxdist250.tsv', sep='\t')
resultTable67x[,'Region']<-factor(resultTable67x[,'Region'])
resultTable67x$Axis <- "x"
resultTable67y <- read.csv(file = 'E6.75 crossing y to t4.5 maxdist250.tsv', sep='\t')
resultTable67y[,'Region']<-factor(resultTable67y[,'Region'])
resultTable67y$Axis <- "y"
resultTable67 <- rbind(resultTable67x,resultTable67y)
resultTable67$Begin.Diff = resultTable67[,c("Begin.Diff")] / MaMuT_scale_factor
resultTable67$End.Diff = resultTable67[,c("End.Diff")] / MaMuT_scale_factor
resultTable67$BE.Diff = resultTable67[,c("End.Diff")] - resultTable67[,c("Begin.Diff")]

#All stages summary dataset
summaryTable <-read_ods(path="all_stages crossing fisher table.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)



#==============
#XY plot of begin and end diff by region
#==============
x_labels <- unique(resultTable65$Region)
new_xy_labels <- c("x" = "Anterior-Posterior", "y" = "Medial-Lateral")
new_13_labels <- c("1" = "", "3" = "")
fills65 <-c("#64ABE3","#FDD8B5")
fills67 <-c("#15B2D1","#DFCE9D")
fills70 <-c("#65CBDA","#F9D199")

resultTable65 %>%
filter( Begin.Diff< (50) ) %>%
  ggplot( aes(x=Begin.Diff, y=End.Diff) ) +

  #geom_tile() +
  #stat_density_2d_filled(geom = "raster",aes(fill = after_stat(density)),contour = FALSE) +
  #stat_density_2d(geom = "polygon", aes(alpha = (..level..) ^ 2, fill = Region)) + 
  #scale_alpha_continuous(range = c(0, 1)) +
    #stat_binhex()+
  #geom_bin2d(bins=100) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(limits=c(0,50), breaks=seq(0, 40, by=10), colours=c("Dark Red","Red","Yellow") )+
  #scale_fill_gradientn(limits=c(0,50), breaks=seq(0, 40, by=10), colours=c("Black","Purple","Yellow") )+
  stat_bin2d(bins=50) +
  scale_fill_viridis() +
  #stat_density2d_filled() +
  #geom_point(size=1, shape=1,alpha=0.1) +
  #geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  #geom_jitter(aes(color=Region), size=4, alpha=0.9) +
  #geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #scale_x_continuous(minor_breaks=c(0)) +
  #scale_y_continuous(minor_breaks=c(0)) +
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  #coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels), labels=c("","","","","","","","","")) +
  #facet_wrap(facets=vars(Region),nrow=1,ncol=2,scales="free_x") +
  stat_cor(aes(label = ..rr.label..), label.y.npc = 0.2, label.x.npc = 0.8, geom="label", vjust=0.5,hjust=0.5, color = "black", size=20) +
  facet_grid( vars(Region),vars(Axis),labeller = labeller(Axis = new_xy_labels,Region=new_13_labels)) +
  #facet_grid( vars(Region))+
  #scale_y_continuous(trans='sqrt') +
  scale_y_continuous(limits=c(-105,155)) +
  scale_x_continuous(limits=c(-1,51)) +
  #scale_y_continuous(limits=c(20,120)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
  #geom_hline(yintercept = 0, color = "#444444") +
  #geom_vline(xintercept = 0, color = "#444444") +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(0, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill=c("Light Blue","White"))
  ) +  
  xlab("")


resultTable67 %>%
  filter( Begin.Diff< (50) ) %>%
  ggplot( aes(x=Begin.Diff, y=End.Diff) ) +
  #geom_tile() +
  #stat_density_2d_filled(geom = "raster",aes(fill = after_stat(density)),contour = FALSE) +
  #stat_density_2d(geom = "polygon", aes(alpha = (..level..) ^ 2, fill = Region)) + 
  #scale_alpha_continuous(range = c(0, 1)) +
  #stat_binhex()+
  #geom_bin2d(bins=100) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(limits=c(0,50), breaks=seq(0, 40, by=10), colours=c("Dark Red","Red","Yellow") )+
  #scale_fill_gradientn(limits=c(0,50), breaks=seq(0, 40, by=10), colours=c("Black","Purple","Yellow") )+
  stat_bin2d(bins=50) +
  scale_fill_viridis() +
  #stat_density2d_filled() +
  #geom_point(size=1, shape=1,alpha=0.1) +
  #geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  #geom_jitter(aes(color=Region), size=4, alpha=0.9) +
  #geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #scale_x_continuous(minor_breaks=c(0)) +
  #scale_y_continuous(minor_breaks=c(0)) +
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  #coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels), labels=c("","","","","","","","","")) +
  #facet_wrap(facets=vars(Region),nrow=1,ncol=2,scales="free_x") +
  stat_cor(aes(label = ..rr.label..), label.y.npc = 0.2, label.x.npc = 0.8, geom="label", vjust=0.5,hjust=0.5, color = "black", size=20) +
  facet_grid( vars(Region),vars(Axis),labeller = labeller(Axis = new_xy_labels,Region=new_13_labels)) +
  #facet_grid( vars(Region))+
  #scale_y_continuous(trans='sqrt') +
  scale_y_continuous(limits=c(-105,155)) +
  scale_x_continuous(limits=c(-1,51)) +
  #scale_y_continuous(limits=c(20,120)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
  #geom_hline(yintercept = 0, color = "#444444") +
  #geom_vline(xintercept = 0, color = "#444444") +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.x = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(0, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill=c("Light Blue","White"))
  ) +  
  xlab("")

#==========
#Density plots
#==========
resultTable67 %>%
filter( Begin.Diff>10 & Begin.Diff<50 ) %>%
  ggplot( aes(x=BE.Diff, fill=Region, color=Region) ) +
  
  #geom_tile() +
  #stat_density_2d_filled(geom = "raster",aes(fill = after_stat(density)),contour = FALSE) +
  #stat_density_2d(geom = "polygon", aes(alpha = (..level..) ^ 2, fill = Region)) + 
  #scale_alpha_continuous(range = c(0, 1)) +
  #stat_binhex()+
  #geom_bin2d(bins=100) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(limits=c(0,50), breaks=seq(0, 40, by=10), colours=c("Dark Red","Red","Yellow") )+
  #scale_fill_gradientn(limits=c(0,50), breaks=seq(0, 40, by=10), colours=c("Black","Purple","Yellow") )+
  #stat_bin2d(bins=50) +
  #scale_fill_viridis() +
  #stat_density2d_filled() +
  #geom_point(size=1, shape=1,alpha=0.1) +
  #geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  #geom_jitter(aes(color=Region), size=4, alpha=0.9) +
  geom_density(adjust=1.5, alpha=0.8) +
  scale_fill_manual(values=fills67) +
  scale_color_manual(values=fills67) +
  #geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #scale_x_continuous(minor_breaks=c(0)) +
  #scale_y_continuous(minor_breaks=c(0)) +
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  #coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels), labels=c("","","","","","","","","")) +
  #facet_wrap(facets=vars(Region),nrow=1,ncol=2,scales="free_x") +
  #stat_cor(aes(label = ..r.label..), label.y.npc = 0.2, label.x.npc = 0.8, geom="label", vjust=0.5,hjust=0.5, color = "black", size=20) +
  #facet_grid( vars(Region),vars(Axis),labeller = labeller(Axis = new_xy_labels,Region=new_13_labels)) +
  facet_wrap( vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #scale_y_continuous(trans='sqrt') +
  #scale_y_continuous(limits=c(-5,5)) +
  scale_x_continuous(limits=c(-70,70)) +
  #scale_y_continuous(limits=c(20,120)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
  #geom_hline(yintercept = 0, color = "#444444") +
  #geom_vline(xintercept = 0, color = "#444444") +
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
    axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White")
  ) +  
  xlab("")

#Density plots
resultTable65 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 ) %>%
  ggplot( aes(x=BE.Diff, fill=Region, color=Region) ) +
  
  #geom_tile() +
  #stat_density_2d_filled(geom = "raster",aes(fill = after_stat(density)),contour = FALSE) +
  #stat_density_2d(geom = "polygon", aes(alpha = (..level..) ^ 2, fill = Region)) + 
  #scale_alpha_continuous(range = c(0, 1)) +
  #stat_binhex()+
  #geom_bin2d(bins=100) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(limits=c(0,50), breaks=seq(0, 40, by=10), colours=c("Dark Red","Red","Yellow") )+
  #scale_fill_gradientn(limits=c(0,50), breaks=seq(0, 40, by=10), colours=c("Black","Purple","Yellow") )+
  #stat_bin2d(bins=50) +
#scale_fill_viridis() +
#stat_density2d_filled() +
#geom_point(size=1, shape=1,alpha=0.1) +
#geom_violin(trim=TRUE,aes(fill = Region),color="black",size=1.5) +
#scale_fill_brewer(palette="Paired") + 
#scale_color_brewer(palette="Paired") +   
#geom_jitter(aes(color=Region), size=4, alpha=0.9) +
geom_density(adjust=1.5, alpha=0.7) +
  scale_fill_manual(values=fills65) +
  scale_color_manual(values=fills65) +
  #geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #scale_x_continuous(minor_breaks=c(0)) +
  #scale_y_continuous(minor_breaks=c(0)) +
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  #coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels), labels=c("","","","","","","","","")) +
  #facet_wrap(facets=vars(Region),nrow=1,ncol=2,scales="free_x") +
  #stat_cor(aes(label = ..r.label..), label.y.npc = 0.2, label.x.npc = 0.8, geom="label", vjust=0.5,hjust=0.5, color = "black", size=20) +
  #facet_grid( vars(Region),vars(Axis),labeller = labeller(Axis = new_xy_labels,Region=new_13_labels)) +
  facet_wrap( vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #scale_y_continuous(trans='sqrt') +
  #scale_y_continuous(limits=c(-5,5)) +
  scale_x_continuous(limits=c(-70,70)) +
  #scale_y_continuous(limits=c(20,120)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
  #geom_hline(yintercept = 0, color = "#444444") +
  #geom_vline(xintercept = 0, color = "#444444") +
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
    axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White")
  ) +  
  xlab("")

#How many data points are there in the density plot
resultTable65 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 ) %>%
  nrow()
resultTable67 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 ) %>%
  nrow()



#Fisher testing
countA <- resultTable65 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 & Region=="1" & BE.Diff >= 0) %>%
  nrow()
countB <- resultTable65 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 & Region=="1" & BE.Diff < 0) %>%
  nrow()
countC <- resultTable65 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 & Region=="3" & BE.Diff >= 0) %>%
  nrow()
countD <- resultTable65 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 & Region=="3" & BE.Diff < 0) %>%
  nrow()
countA[1]
countB[1]
countC[1]
countD[1]
fisher.test(matrix(unlist(c(countA[1], countB[1], countC[1], countD[1])), nrow = 2))

countA <- resultTable67 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 & Region=="1" & BE.Diff >= 0) %>%
  nrow()
countB <- resultTable67 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 & Region=="1" & BE.Diff < 0) %>%
  nrow()
countC <- resultTable67 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 & Region=="3" & BE.Diff >= 0) %>%
  nrow()
countD <- resultTable67 %>%
  filter( Begin.Diff>10 & Begin.Diff<50 & Region=="3" & BE.Diff < 0) %>%
  nrow()
countA[1]
countB[1]
countC[1]
countD[1]
fisher.test(matrix(unlist(c(countA[1], countB[1], countC[1], countD[1])), nrow = 2))


#==============
#summary table plots
#==============
x_labels <- c( "E6.5", "E6.75", "E7.0")
sand <- c(fills65[2],fills67[2],fills70[2])
water <- c(fills65[1],fills67[1],fills70[1])
summaryTable %>%
  ggplot( aes(x=Stage, y=AvgNumCrosses, fill=Stage, color=Stage, alpha=0.4)) +
  geom_boxplot( outlier.shape = NA, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=10, alpha=0.8) +
  coord_flip() + 
  facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  scale_y_continuous(trans="log") +
  scale_x_discrete(labels= x_labels, limits = rev) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=water) +
  scale_fill_manual(values=sand) +
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
    axis.text.y = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
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

