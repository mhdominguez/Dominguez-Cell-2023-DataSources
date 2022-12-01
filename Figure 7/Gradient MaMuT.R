#install.packages('readODS')
#font_import()
library(ggplot2)
library(readODS)
library(tidyverse)
library(viridis)
library(dplyr)
library(extrafont)
library(ggridges)
library(patchwork)
library(circular)
loadfonts()

library("RColorBrewer")
display.brewer.all()
brewer.pal(n = 9, name = "Paired")

#basically, if going 0->30 versus 60->90, pretty much have to use den12 to get gradient inversion
#if going 0->60 and 60->120, can use either den5 or den12
#have to do X+Y always, never 2X+Y

quant_breaks = 9
quant_breaks_radial_histogram = 3
#XY_remove_cutoff = 100 #remove extreme 10% of X extent and Y extent, to try to democratize density measurements


#===================================================================================
#
#
# E6.75 KO vs control
#
#
#===================================================================================

################
# t0 - t120
################

#open E6.75 KO dataset
resultTable_KO <- read.csv(file = "E6.75 KO data.tsv", sep="\t" )

concentric_ring_gradient_KO <- function() {
  #set up coordinate system and concentric ring gradient
  X_start = min(resultTable_KO$End.X)
  Y_start = min(resultTable_KO$End.Y)
  resultTable_KO$End.X <- resultTable_KO$End.X - X_start
  resultTable_KO$End.Y <- resultTable_KO$End.Y - Y_start
  resultTable_KO$Begin.X <- resultTable_KO$Begin.X - X_start
  resultTable_KO$Begin.Y <- resultTable_KO$Begin.Y - Y_start
  resultTable_KO$End.XYComp <- sqrt( (resultTable_KO$End.Y)^2 + (resultTable_KO$End.X)^2 )
  resultTable_KO$End.XYComp.Cut_values <- cut_number(resultTable_KO$End.XYComp,n=quant_breaks_radial_histogram)
  resultTable_KO$End.XYComp.Cut_factor <- as_factor(as.numeric(resultTable_KO$End.XYComp.Cut_values))
  quant_KO <- quantile(resultTable_KO$End.XYComp, probs = seq(0, 1, 1/quant_breaks) )
  
  resultTable_KO$Condition <- "KO"
  
  #remove extremes along either axis to mitigate density dropoffs
  #X_start = 0
  #X_stop = max(resultTable_KO$End.X)
  #X_range = (X_stop-X_start) / XY_remove_cutoff
  #X_bottom_KO = X_start + X_range
  #X_top_KO = X_stop - X_range
  #Y_start = 0
  #Y_stop = max(resultTable_KO$End.Y)
  #Y_range = (Y_stop-Y_start) / XY_remove_cutoff
  #Y_bottom_KO = Y_start + Y_range
  #Y_top_KO = Y_stop - Y_range
  
  #add direction data, XY plane only
  resultTable_KO$Track.Vector.Angle <- 180 * (atan2((resultTable_KO$End.Y-resultTable_KO$Begin.Y),(resultTable_KO$End.X-resultTable_KO$Begin.X))) / pi + 180
  
  resultTable_KO <<- resultTable_KO
  quant_KO <<- quant_KO
}

linear_XY_gradient_KO <- function() {
  #Ymed = median(resultTable_KO$End.Y)
  #Xmed = median(resultTable_KO$End.X)
  #resultTable_KO$End.X <- resultTable_KO$End.X - Xmed
  #resultTable_KO$End.Y <- resultTable_KO$End.Y - Ymed
  #resultTable_KO$Begin.X <- resultTable_KO$Begin.X - Xmed
  #resultTable_KO$Begin.Y <- resultTable_KO$Begin.Y - Ymed
  ##resultTable_KO$Begin.X <- - resultTable_KO$Begin.X
  ##resultTable_KO$End.X <- - resultTable_KO$End.X
  resultTable_KO$End.XYComp <- (1* resultTable_KO$End.Y) + resultTable_KO$End.X
  resultTable_KO$End.XYComp.Cut_values <- cut_number(resultTable_KO$End.XYComp,n=quant_breaks_radial_histogram)
  resultTable_KO$End.XYComp.Cut_factor <- as_factor(as.numeric(resultTable_KO$End.XYComp.Cut_values))
  resultTable_KO$Condition <- "KO"
  quant_KO <- quantile(resultTable_KO$End.XYComp, probs = seq(0, 1, 1/quant_breaks) )
  
  #remove extremes along either axis to mitigate density dropoffs
  #X_start = min(resultTable_KO$End.X)
  #X_stop = max(resultTable_KO$End.X)
  #X_range = (X_stop-X_start) / XY_remove_cutoff
  #X_bottom_KO = X_start + X_range
  #X_top_KO = X_stop - X_range
  #Y_start = min(resultTable_KO$End.Y)
  #Y_stop = max(resultTable_KO$End.Y)
  #Y_range = (Y_stop-Y_start) / XY_remove_cutoff
  #Y_bottom_KO = Y_start + Y_range
  #Y_top_KO = Y_stop - Y_range
  
  #add direction data, XY plane only
  resultTable_KO$Track.Vector.Angle <- 180 * (atan2((resultTable_KO$End.Y-resultTable_KO$Begin.Y),(resultTable_KO$End.X-resultTable_KO$Begin.X))) / pi + 180

  resultTable_KO <<- resultTable_KO
  quant_KO <<- quant_KO
}

linear_XY_gradient_KO() 
#concentric_ring_gradient_KO()

#open E6.75 ctrl dataset
resultTable_ctrl <- read.csv(file = "E6.75 control data.tsv", sep="\t" )

concentric_ring_gradient_ctrl <- function() {
  #reverse x axis since embryo flipped
  Xmed = median(resultTable_ctrl$End.X)
  resultTable_ctrl$End.X <- resultTable_ctrl$End.X - Xmed
  resultTable_ctrl$Begin.X <- resultTable_ctrl$Begin.X - Xmed
  resultTable_ctrl$Begin.X <- - resultTable_ctrl$Begin.X
  resultTable_ctrl$End.X <- - resultTable_ctrl$End.X
  
  #set up coordinate system and concentric ring gradient
  X_start = min(resultTable_ctrl$End.X)
  Y_start = min(resultTable_ctrl$End.Y)
  resultTable_ctrl$End.X <- resultTable_ctrl$End.X - X_start
  resultTable_ctrl$End.Y <- resultTable_ctrl$End.Y - Y_start
  resultTable_ctrl$Begin.X <- resultTable_ctrl$Begin.X - X_start
  resultTable_ctrl$Begin.Y <- resultTable_ctrl$Begin.Y - Y_start
  resultTable_ctrl$End.XYComp <- sqrt( (resultTable_ctrl$End.Y)^2 + (resultTable_ctrl$End.X)^2 )
  quant_ctrl <- quantile(resultTable_ctrl$End.XYComp, probs = seq(0, 1, 1/quant_breaks) )
  resultTable_ctrl$End.XYComp.Cut_values <- cut_number(resultTable_ctrl$End.XYComp,n=quant_breaks_radial_histogram)
  resultTable_ctrl$End.XYComp.Cut_factor <- as_factor(as.numeric(resultTable_ctrl$End.XYComp.Cut_values))
  resultTable_ctrl$Condition <- "Het"
  
  #remove extremes along either axis to mitigate density dropoffs
  #X_start = 0
  #X_stop = max(resultTable_ctrl$End.X)
  #X_range = (X_stop-X_start) / XY_remove_cutoff
  #X_bottom_ctrl = X_start + X_range
  #X_top_ctrl = X_stop - X_range
  #Y_start = 0
  #Y_stop = max(resultTable_ctrl$End.Y)
  #Y_range = (Y_stop-Y_start) / XY_remove_cutoff
  #Y_bottom_ctrl = Y_start + Y_range
  #Y_top_ctrl = Y_stop - Y_range
  
  #add direction data, XY plane only
  resultTable_ctrl$Track.Vector.Angle <- 180 * (atan2((resultTable_ctrl$End.Y-resultTable_ctrl$Begin.Y),(resultTable_ctrl$End.X-resultTable_ctrl$Begin.X))) / pi + 180
  
  resultTable_ctrl <<- resultTable_ctrl
  quant_ctrl <<- quant_ctrl
}

linear_XY_gradient_ctrl <- function() {
  Ymed = median(resultTable_ctrl$End.Y)
  #Xmed = median(resultTable_ctrl$End.X)
  #resultTable_ctrl$End.X <- resultTable_ctrl$End.X - Xmed
  #resultTable_ctrl$End.Y <- resultTable_ctrl$End.Y - Ymed
  #resultTable_ctrl$Begin.X <- resultTable_ctrl$Begin.X - Xmed
  #resultTable_ctrl$Begind.Y <- resultTable_ctrl$Begin.Y - Ymed
  resultTable_ctrl$Begin.X <- - resultTable_ctrl$Begin.X
  resultTable_ctrl$End.X <- - resultTable_ctrl$End.X
  resultTable_ctrl$End.XYComp <- (1 * resultTable_ctrl$End.Y ) + resultTable_ctrl$End.X
  resultTable_ctrl$End.XYComp.Cut_values <- cut_number(resultTable_ctrl$End.XYComp,n=quant_breaks_radial_histogram)
  resultTable_ctrl$End.XYComp.Cut_factor <- as_factor(as.numeric(resultTable_ctrl$End.XYComp.Cut_values))
  resultTable_ctrl$Condition <- "Het"
  quant_ctrl <- quantile(resultTable_ctrl$End.XYComp, probs = seq(0, 1, 1/quant_breaks) )
  
  resultTable_ctrl$Condition <- "Het"
  
  #remove extremes along either axis to mitigate density dropoffs
  #X_start = min(resultTable_ctrl$End.X)
  #X_stop = max(resultTable_ctrl$End.X)
  #X_range = (X_stop-X_start) / XY_remove_cutoff
  #X_bottom_ctrl = X_start + X_range
  #X_top_ctrl = X_stop - X_range
  #Y_start = min(resultTable_ctrl$End.Y)
  #Y_stop = max(resultTable_ctrl$End.Y)
  #Y_range = (Y_stop-Y_start) / XY_remove_cutoff
  #Y_bottom_ctrl = Y_start + Y_range
  #Y_top_ctrl = Y_stop - Y_range
  
  #add direction data, XY plane only
  resultTable_ctrl$Track.Vector.Angle <- 180 * (atan2((resultTable_ctrl$End.Y-resultTable_ctrl$Begin.Y),(resultTable_ctrl$End.X-resultTable_ctrl$Begin.X))) / pi + 180
  
  resultTable_ctrl <<- resultTable_ctrl
  quant_ctrl <<- quant_ctrl
}

linear_XY_gradient_ctrl()
#concentric_ring_gradient_ctrl()

#combine the datasets
resultTable_combined <- rbind(resultTable_ctrl,resultTable_KO )
resultTable_combined$Condition <- as_factor(resultTable_combined$Condition)

################
# Generate t0 - t120 density plots
################
p1 <- resultTable_ctrl %>%
#  filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
 # filter( End.X > X_bottom_ctrl & End.X < X_top_ctrl & End.Y > Y_bottom_ctrl & End.Y < Y_top_ctrl ) %>%
  filter( End.XYComp > quant_ctrl[2] & End.XYComp < quant_ctrl[quant_breaks] ) %>%
  ggplot( aes(x=End.XYComp, y=Avg.Density..cells.per.12.radii., alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  #scale_fill_brewer(palette="RdYlBu") + 
  #scale_color_brewer(palette="RdYlBu") +  
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
#  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
 geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_ctrl,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) + #,limits=c(500,2000)) +
#  scale_x_binned(n.breaks=7) + #,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) +
  scale_y_continuous(limits=c(1,150),trans='sqrt') +
  #scale_x_binned() +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
 # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
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

p2 <- resultTable_KO %>%
 # filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
 # filter( End.X > X_bottom_KO & End.X < X_top_KO & End.Y > Y_bottom_KO, End.Y < Y_top_KO ) %>%
  #filter( Begin.T > 60 ) %>%
  filter( End.XYComp > quant_KO[2] & End.XYComp < quant_KO[quant_breaks] ) %>%
  ggplot( aes(x=End.XYComp, y=Avg.Density..cells.per.12.radii., alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_KO,limits=c(quant_KO[2],quant_KO[quant_breaks])) + #,limits=c(500,2000)) +
 #scale_x_binned(n.breaks=7) + #,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) +
  scale_y_continuous(limits=c(1,150),trans='sqrt') +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
 # scale_y_continuous(limits=c(1,21),trans='sqrt') +
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

p1+p2

test_3 <- filter(resultTable_ctrl, End.XYComp > quant_ctrl[2] & End.XYComp < quant_ctrl[3] )
test_4 <- filter(resultTable_ctrl, End.XYComp > quant_ctrl[quant_breaks-1] & End.XYComp < quant_ctrl[quant_breaks] )
t.test(test_3$Avg.Density..cells.per.12.radii.,test_4$Avg.Density..cells.per.12.radii.)

test_1 <- filter(resultTable_KO, End.XYComp > quant_KO[2] & End.XYComp < quant_KO[3] )
test_2 <- filter(resultTable_KO, End.XYComp > quant_KO[quant_breaks-1] & End.XYComp < quant_KO[quant_breaks] )
t.test(test_1$Avg.Density..cells.per.12.radii.,test_2$Avg.Density..cells.per.12.radii.)




################
# Generate t0 - t120 motility plots
################
p1 <- resultTable_ctrl %>%
  filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_ctrl & End.X < X_top_ctrl & End.Y > Y_bottom_ctrl & End.Y < Y_top_ctrl ) %>%
  # filter( Avg.Sliding.Velocity..micron.hr. > 20 ) %>%
  filter( End.XYComp > quant_ctrl[2] & End.XYComp < quant_ctrl[quant_breaks] ) %>%  
  ggplot( aes(x=End.XYComp, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  #scale_fill_brewer(palette="RdYlBu") + 
  #scale_color_brewer(palette="RdYlBu") +  
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_ctrl,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) + #,limits=c(500,2000)) +
  #  scale_x_binned(n.breaks=7) + #,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) +
  scale_y_continuous(trans='sqrt',limits=c(11,101)) +
  #scale_x_binned() +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
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

p2 <- resultTable_KO %>%
  filter( End.XYComp > quant_KO[2] & End.XYComp < quant_KO[quant_breaks] ) %>%
  filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_KO & End.X < X_top_KO & End.Y > Y_bottom_KO, End.Y < Y_top_KO ) %>%
  #filter( Begin.T > 60 ) %>%
  ggplot( aes(x=End.XYComp, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_KO,limits=c(quant_KO[2],quant_KO[quant_breaks])) + #,limits=c(500,2000)) +
  #scale_x_binned(n.breaks=7) + #,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
  scale_y_continuous(limits=c(11,101),trans='sqrt') +
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

p1+p2

test_3 <- filter(resultTable_ctrl, End.XYComp > quant_ctrl[2] & End.XYComp < quant_ctrl[3] & Avg.Sliding.Velocity..micron.hr. < 100 )
test_4 <- filter(resultTable_ctrl, End.XYComp > quant_ctrl[quant_breaks-1] & End.XYComp < quant_ctrl[quant_breaks] & Avg.Sliding.Velocity..micron.hr. < 100 )
t.test(test_3$Avg.Sliding.Velocity..micron.hr.,test_4$Avg.Sliding.Velocity..micron.hr.)

test_1 <- filter(resultTable_KO, End.XYComp > quant_KO[2] & End.XYComp < quant_KO[3] & Avg.Sliding.Velocity..micron.hr. < 100 )
test_2 <- filter(resultTable_KO, End.XYComp > quant_KO[quant_breaks-1] & End.XYComp < quant_KO[quant_breaks] & Avg.Sliding.Velocity..micron.hr. < 100 )
t.test(test_1$Avg.Sliding.Velocity..micron.hr.,test_2$Avg.Sliding.Velocity..micron.hr.)

t.test(test_4$Avg.Sliding.Velocity..micron.hr.,test_2$Avg.Sliding.Velocity..micron.hr.)
t.test(test_1$Avg.Sliding.Velocity..micron.hr.,test_3$Avg.Sliding.Velocity..micron.hr.)

t.test(resultTable_ctrl$Avg.Sliding.Velocity..micron.hr.,resultTable_KO$Avg.Sliding.Velocity..micron.hr.)



################
# Generate t0 - t120 birthday plot
################

p1 <- resultTable_ctrl %>%
  filter( End.XYComp > quant_ctrl[2] & End.XYComp < quant_ctrl[quant_breaks] ) %>%
  #filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_ctrl & End.X < X_top_ctrl & End.Y > Y_bottom_ctrl & End.Y < Y_top_ctrl ) %>%
  # filter( Avg.Sliding.Velocity..micron.hr. > 20 ) %>%
  ggplot( aes(x=End.XYComp, y=Begin.T, alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  #scale_fill_brewer(palette="RdYlBu") + 
  #scale_color_brewer(palette="RdYlBu") +  
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_ctrl,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) + #,limits=c(500,2000)) +
  #  scale_x_binned(n.breaks=7) + #,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) +
 # scale_y_continuous(limits=c(1,21),trans='sqrt') +
  #scale_x_binned() +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
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

p2 <- resultTable_KO %>%
  filter( End.XYComp > quant_KO[2] & End.XYComp < quant_KO[quant_breaks] ) %>%
  # filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_KO & End.X < X_top_KO & End.Y > Y_bottom_KO, End.Y < Y_top_KO ) %>%
  #filter( Begin.T > 60 ) %>%
  ggplot( aes(x=End.XYComp, y=Begin.T, alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_KO,limits=c(quant_KO[2],quant_KO[quant_breaks])) + #,limits=c(500,2000)) +
  #scale_x_binned(n.breaks=7) + #,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
  #scale_y_continuous(limits=c(1,21),trans='sqrt') +
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

p1+p2

test_3 <- filter(resultTable_ctrl, End.XYComp > quant_ctrl[2] & End.XYComp < quant_ctrl[3] )
test_4 <- filter(resultTable_ctrl, End.XYComp > quant_ctrl[quant_breaks-1] & End.XYComp < quant_ctrl[quant_breaks] )
wilcox.test(test_3$Begin.T,test_4$Begin.T)

test_1 <- filter(resultTable_KO, End.XYComp > quant_KO[2] & End.XYComp < quant_KO[3] )
test_2 <- filter(resultTable_KO, End.XYComp > quant_KO[quant_breaks-1] & End.XYComp < quant_KO[quant_breaks] )
wilcox.test(test_1$Begin.T,test_2$Begin.T)



#==============
#radial histogram of angles for trajectories of tracks
#==============
resultTable_combined %>%
  #filter( End.X > X_bottom_KO & End.X < X_top_KO & End.Y > Y_bottom_KO, End.Y < Y_top_KO ) %>%
  
  ggplot( aes(x=End.XYComp, fill=Condition, color=Condition, alpha=0.8)) +
  #ggplot( aes(x=End.XYComp, y=Track.Vector.Angle, alpha=0.6, color="Black")) +
  geom_density(aes(x=Track.Vector.Angle,fill=Condition),alpha=0.5,binwidth=10,size=2 ) + #,position = "stack") + 
  facet_wrap(~End.XYComp.Cut_factor, ncol = quant_breaks_radial_histogram) +
  #scale_x_continuous(breaks=seq(15,360,30),limits=c(0,360)) + 
  #scale_x_continuous(expand=c(0,-5),limits=c(0,360),breaks=seq(5,365,30)) +
  scale_x_continuous(limits=c(0,360),breaks=seq(5,365,30)) +
  scale_y_continuous(trans='sqrt') +
  coord_polar(theta="x", start=3*pi/2-0.0872665, direction=-1) + 
  #theme_bw() +
  #geom_boxplot( outlier.shape = NA, lwd=2 ) +
  #geom_crossbar(stat = "summary", fun=mean,color="black",size=2.5) +
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
  scale_color_manual(values=c("#000000", "#800000" )) + #, "#808080", "#808080", "#808080", "#808080")) +
  #scale_color_manual(values=water) +
  scale_fill_manual(values=c("#000000", "#F98A81" )) + 
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
    panel.grid.major.x = element_line(colour="#707070",size=0.6),
    panel.grid.major.y = element_line(colour="#707070",size=1.2),
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

test_1 <- filter(resultTable_combined, str_detect(Condition,"Het"), End.XYComp.Cut_factor == 1 )
test_2 <- filter(resultTable_combined, str_detect(Condition,"KO"),End.XYComp.Cut_factor == 1 )
test_1_circular <- as.circular(test_1$Track.Vector.Angle,units="degrees")
test_2_circular <- as.circular(test_2$Track.Vector.Angle,units="degrees")

watson.two.test(test_1_circular,test_2_circular) #,alpha=0.05)


test_1 <- filter(resultTable_combined, str_detect(Condition,"Het"), End.XYComp.Cut_factor == 2 )
test_2 <- filter(resultTable_combined, str_detect(Condition,"KO"),End.XYComp.Cut_factor == 2 )
test_1_circular <- as.circular(test_1$Track.Vector.Angle,units="degrees")
test_2_circular <- as.circular(test_2$Track.Vector.Angle,units="degrees")

watson.two.test(test_1_circular,test_2_circular) #,alpha=0.05)


test_1 <- filter(resultTable_combined, str_detect(Condition,"Het"), End.XYComp.Cut_factor == quant_breaks_radial_histogram )
test_2 <- filter(resultTable_combined, str_detect(Condition,"KO"),End.XYComp.Cut_factor == quant_breaks_radial_histogram )
test_1_circular <- as.circular(test_1$Track.Vector.Angle,units="degrees")
test_2_circular <- as.circular(test_2$Track.Vector.Angle,units="degrees")

watson.two.test(test_1_circular,test_2_circular) #,alpha=0.05)


test_1 <- filter(resultTable_combined, str_detect(Condition,"Het"), End.XYComp.Cut_factor == quant_breaks_radial_histogram-1 )
test_2 <- filter(resultTable_combined, str_detect(Condition,"KO"),End.XYComp.Cut_factor == quant_breaks_radial_histogram-1 )
test_1_circular <- as.circular(test_1$Track.Vector.Angle,units="degrees")
test_2_circular <- as.circular(test_2$Track.Vector.Angle,units="degrees")

watson.two.test(test_1_circular,test_2_circular) #,alpha=0.05)



#===================================================================================
#
#
# control E6.5 vs. E6.75 vs. E7.0
#
#
#===================================================================================

resultTable_E6.75 <- resultTable_ctrl
quant_E6.75 <- quant_ctrl

#open E6.5 ctrl dataset
resultTable_ctrl <- read.csv(file = "E6.5 control data.tsv", sep="\t" )

linear_XY_gradient_ctrl()

resultTable_E6.5 <- resultTable_ctrl
quant_E6.5 <- quant_ctrl

#open E7.0 ctrl dataset
resultTable_ctrl <- read.csv(file = "E7.0 control data.tsv", sep="\t" )

linear_XY_gradient_ctrl()

resultTable_E7.0 <- resultTable_ctrl
quant_E7.0 <- quant_ctrl


################
# Generate t0 - t120 birthday plot
################

p1 <- resultTable_E6.5 %>%
  filter( End.XYComp > quant_E6.5[2] & End.XYComp < quant_E6.5[quant_breaks] ) %>%
  #filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_ctrl & End.X < X_top_ctrl & End.Y > Y_bottom_ctrl & End.Y < Y_top_ctrl ) %>%
  # filter( Avg.Sliding.Velocity..micron.hr. > 20 ) %>%
  ggplot( aes(x=End.XYComp, y=Begin.T, alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  #scale_fill_brewer(palette="RdYlBu") + 
  #scale_color_brewer(palette="RdYlBu") +  
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_E6.5,limits=c(quant_E6.5[2],quant_E6.5[quant_breaks])) + #,limits=c(500,2000)) +
  #  scale_x_binned(n.breaks=7) + #,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) +
  # scale_y_continuous(limits=c(1,21),trans='sqrt') +
  #scale_x_binned() +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
 # scale_y_continuous(trans='sqrt')+
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

p2 <- resultTable_E6.75 %>%
  filter( End.XYComp > quant_E6.75[2] & End.XYComp < quant_E6.75[quant_breaks] ) %>%
  # filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_KO & End.X < X_top_KO & End.Y > Y_bottom_KO, End.Y < Y_top_KO ) %>%
  #filter( Begin.T > 60 ) %>%
  ggplot( aes(x=End.XYComp, y=Begin.T, alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_E6.75,limits=c(quant_E6.75[2],quant_E6.75[quant_breaks])) + #,limits=c(500,2000)) +
  #scale_x_binned(n.breaks=7) + #,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
  #scale_y_continuous(limits=c(1,21),trans='sqrt') +
  #scale_y_continuous(limits=c(20,120)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
  #scale_y_continuous(trans='sqrt')+
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
p3 <- resultTable_E7.0 %>%
  filter( End.XYComp > quant_E7.0[2] & End.XYComp < quant_E7.0[quant_breaks] ) %>%
  # filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_KO & End.X < X_top_KO & End.Y > Y_bottom_KO, End.Y < Y_top_KO ) %>%
  #filter( Begin.T > 60 ) %>%
  ggplot( aes(x=End.XYComp, y=Begin.T, alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_E7.0,limits=c(quant_E7.0[2],quant_E7.0[quant_breaks])) + #,limits=c(500,2000)) +
  #scale_x_binned(n.breaks=7) + #,limits=c(quant_ctrl[2],quant_ctrl[quant_breaks])) +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
  #scale_y_continuous(limits=c(1,21),trans='sqrt') +
  #scale_y_continuous(limits=c(20,120)) +
  #scale_fill_grey(start=0.0,end=0.0) +
  #scale_fill_manual(values=c("red", "dark red")) +
 # scale_y_continuous(trans='sqrt')+
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
p1+p2+p3

test_3 <- filter(resultTable_E6.5, End.XYComp > quant_E6.5[2] & End.XYComp < quant_E6.5[3] )
test_4 <- filter(resultTable_E6.5, End.XYComp > quant_E6.5[quant_breaks-1] & End.XYComp < quant_E6.5[quant_breaks] )
t.test(test_3$Begin.T,test_4$Begin.T)

test_1 <- filter(resultTable_E6.75, End.XYComp > quant_E6.75[2] & End.XYComp < quant_E6.75[3] )
test_2 <- filter(resultTable_E6.75, End.XYComp > quant_E6.75[quant_breaks-1] & End.XYComp < quant_E6.75[quant_breaks] )
t.test(test_1$Begin.T,test_2$Begin.T)

test_5 <- filter(resultTable_E7.0, End.XYComp > quant_E7.0[2] & End.XYComp < quant_E7.0[3] )
test_6 <- filter(resultTable_E7.0, End.XYComp > quant_E7.0[quant_breaks-1] & End.XYComp < quant_E7.0[quant_breaks] )
t.test(test_5$Begin.T,test_6$Begin.T)


################
# Generate t0 - t120 motility plot
################
p1 <- resultTable_E6.5 %>%
  filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_E6.5 & End.X < X_top_E6.5 & End.Y > Y_bottom_E6.5 & End.Y < Y_top_E6.5 ) %>%
  # filter( Avg.Sliding.Velocity..micron.hr. > 20 ) %>%
  filter( End.XYComp > quant_E6.5[2] & End.XYComp < quant_E6.5[quant_breaks] ) %>%  
  ggplot( aes(x=End.XYComp, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  #scale_fill_brewer(palette="RdYlBu") + 
  #scale_color_brewer(palette="RdYlBu") +  
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_E6.5,limits=c(quant_E6.5[2],quant_E6.5[quant_breaks])) + #,limits=c(500,2000)) +
  #  scale_x_binned(n.breaks=7) + #,limits=c(quant_E6.5[2],quant_E6.5[quant_breaks])) +
  scale_y_continuous(trans='sqrt',limits=c(11,101)) +
  #scale_x_binned() +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
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

p2 <- resultTable_E6.75 %>%
  filter( End.XYComp > quant_E6.75[2] & End.XYComp < quant_E6.75[quant_breaks] ) %>%
  filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_E6.75 & End.X < X_top_E6.75 & End.Y > Y_bottom_E6.75, End.Y < Y_top_E6.75 ) %>%
  #filter( Begin.T > 60 ) %>%
  ggplot( aes(x=End.XYComp, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_E6.75,limits=c(quant_E6.75[2],quant_E6.75[quant_breaks])) + #,limits=c(500,2000)) +
  #scale_x_binned(n.breaks=7) + #,limits=c(quant_E6.5[2],quant_E6.5[quant_breaks])) +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
  scale_y_continuous(limits=c(11,101),trans='sqrt') +
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

p3 <- resultTable_E7.0 %>%
  filter( End.XYComp > quant_E7.0[2] & End.XYComp < quant_E7.0[quant_breaks] ) %>%
  filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_E7.0 & End.X < X_top_E7.0 & End.Y > Y_bottom_E7.0, End.Y < Y_top_E7.0 ) %>%
  #filter( Begin.T > 60 ) %>%
  ggplot( aes(x=End.XYComp, y=Avg.Sliding.Velocity..micron.hr., alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_E7.0,limits=c(quant_E7.0[2],quant_E7.0[quant_breaks])) + #,limits=c(500,2000)) +
  #scale_x_binned(n.breaks=7) + #,limits=c(quant_E6.5[2],quant_E6.5[quant_breaks])) +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
  scale_y_continuous(limits=c(11,101),trans='sqrt') +
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

p1+p2+p3

test_3 <- filter(resultTable_E6.5, End.XYComp > quant_E6.5[2] & End.XYComp < quant_E6.5[3] & Avg.Sliding.Velocity..micron.hr. < 100)
test_4 <- filter(resultTable_E6.5, End.XYComp > quant_E6.5[quant_breaks-1] & End.XYComp < quant_E6.5[quant_breaks] & Avg.Sliding.Velocity..micron.hr. < 100 )
t.test(test_3$Avg.Sliding.Velocity..micron.hr.,test_4$Avg.Sliding.Velocity..micron.hr.)

test_1 <- filter(resultTable_E6.75, End.XYComp > quant_E6.75[2] & End.XYComp < quant_E6.75[3]  & Avg.Sliding.Velocity..micron.hr. < 100)
test_2 <- filter(resultTable_E6.75, End.XYComp > quant_E6.75[quant_breaks-1] & End.XYComp < quant_E6.75[quant_breaks] & Avg.Sliding.Velocity..micron.hr. < 100 )
t.test(test_1$Avg.Sliding.Velocity..micron.hr.,test_2$Avg.Sliding.Velocity..micron.hr.)

test_5 <- filter(resultTable_E7.0, End.XYComp > quant_E7.0[2] & End.XYComp < quant_E7.0[3]  & Avg.Sliding.Velocity..micron.hr. < 100)
test_6 <- filter(resultTable_E7.0, End.XYComp > quant_E7.0[quant_breaks-1] & End.XYComp < quant_E7.0[quant_breaks] & Avg.Sliding.Velocity..micron.hr. < 100 )
t.test(test_5$Avg.Sliding.Velocity..micron.hr.,test_6$Avg.Sliding.Velocity..micron.hr.)



################
# Generate t0 - t120 density plots
################
p1 <- resultTable_E6.5 %>%
  #  filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_E6.5 & End.X < X_top_E6.5 & End.Y > Y_bottom_E6.5 & End.Y < Y_top_E6.5 ) %>%
  filter( End.XYComp > quant_E6.5[2] & End.XYComp < quant_E6.5[quant_breaks] ) %>%
  ggplot( aes(x=End.XYComp, y=Avg.Density..cells.per.12.radii., alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  #scale_fill_brewer(palette="RdYlBu") + 
  #scale_color_brewer(palette="RdYlBu") +  
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_E6.5,limits=c(quant_E6.5[2],quant_E6.5[quant_breaks])) + #,limits=c(500,2000)) +
  #  scale_x_binned(n.breaks=7) + #,limits=c(quant_E6.5[2],quant_E6.5[quant_breaks])) +
  scale_y_continuous(limits=c(1,150),trans='sqrt') +
  #scale_x_binned() +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
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

p2 <- resultTable_E6.75 %>%
  # filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_E6.75 & End.X < X_top_E6.75 & End.Y > Y_bottom_E6.75, End.Y < Y_top_E6.75 ) %>%
  #filter( Begin.T > 60 ) %>%
  filter( End.XYComp > quant_E6.75[2] & End.XYComp < quant_E6.75[quant_breaks] ) %>%
  ggplot( aes(x=End.XYComp, y=Avg.Density..cells.per.12.radii., alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_E6.75,limits=c(quant_E6.75[2],quant_E6.75[quant_breaks])) + #,limits=c(500,2000)) +
  #scale_x_binned(n.breaks=7) + #,limits=c(quant_E6.5[2],quant_E6.5[quant_breaks])) +
  scale_y_continuous(limits=c(1,150),trans='sqrt') +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
  # scale_y_continuous(limits=c(1,21),trans='sqrt') +
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

p3 <- resultTable_E7.0 %>%
  # filter( Avg.Sliding.Velocity..micron.hr. < 100 ) %>%
  # filter( End.X > X_bottom_E7.0 & End.X < X_top_E7.0 & End.Y > Y_bottom_E7.0, End.Y < Y_top_E7.0 ) %>%
  #filter( Begin.T > 60 ) %>%
  filter( End.XYComp > quant_E7.0[2] & End.XYComp < quant_E7.0[quant_breaks] ) %>%
  ggplot( aes(x=End.XYComp, y=Avg.Density..cells.per.12.radii., alpha=0.6, color="Black")) +
  #geom_point(aes(fill = Region),color="black",size=1.5) +
  #geom_bar() +
  #scale_fill_brewer(palette="Paired") + 
  #scale_color_brewer(palette="Paired") +   
  scale_color_gradientn(colours = c("Green","Cyan","Blue")) +
  scale_fill_gradientn(colours = c("Green","Cyan","Blue")) +
  geom_jitter(aes(color=End.XYComp), size=4, alpha=0.6) +
  #  stat_density_2d_filled() + #(aes(color="Black")) +
  #geom_crossbar(aes(color=Region), stat = "summary", fun=mean,size=2.5) +
  geom_crossbar(color="Black", stat = "summary", fun=mean,size=2.5) +
  #scale_x_discrete(labels=c("Upper Left","Upper Middle","Upper Right","Middle Left","Middle Middle","Middle Right","Lower Left","Lower Middle","Lower Right")) +
  ylab("") +  
  #ggtitle(expression(bold(paste(bolditalic("Rosa26"^{"Ai66"}), " by 3d region")))) +
  #theme_dark() +
  #coord_cartesian(ylim=c(0, 120)) +
  scale_x_binned(breaks=quant_E7.0,limits=c(quant_E7.0[2],quant_E7.0[quant_breaks])) + #,limits=c(500,2000)) +
  #scale_x_binned(n.breaks=7) + #,limits=c(quant_E6.5[2],quant_E6.5[quant_breaks])) +
  scale_y_continuous(limits=c(1,150),trans='sqrt') +
  #scale_x_continuous(limits=c(500,2000)) +
  coord_flip() + 
  #scale_x_discrete(limits = rev(x_labels)) + #, labels=c("JCF","aFHF","aSHF")) +
  #scale_y_continuous(limits=c(0,150)) +
  #facet_wrap(facets=vars(Region),nrow=3,ncol=3)+#,scales="free_x") +
  # scale_y_continuous(trans='sqrt',limits=c(0,125)) +
  # scale_y_continuous(limits=c(1,21),trans='sqrt') +
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
p1+p2+p3

test_3 <- filter(resultTable_E6.5, End.XYComp > quant_E6.5[2] & End.XYComp < quant_E6.5[3] )
test_4 <- filter(resultTable_E6.5, End.XYComp > quant_E6.5[quant_breaks-1] & End.XYComp < quant_E6.5[quant_breaks] )
t.test(test_3$Avg.Density..cells.per.12.radii.,test_4$Avg.Density..cells.per.12.radii.)

test_1 <- filter(resultTable_E6.75, End.XYComp > quant_E6.75[2] & End.XYComp < quant_E6.75[3] )
test_2 <- filter(resultTable_E6.75, End.XYComp > quant_E6.75[quant_breaks-1] & End.XYComp < quant_E6.75[quant_breaks] )
t.test(test_1$Avg.Density..cells.per.12.radii.,test_2$Avg.Density..cells.per.12.radii.)

test_5 <- filter(resultTable_E7.0, End.XYComp > quant_E7.0[2] & End.XYComp < quant_E7.0[3] )
test_6 <- filter(resultTable_E7.0, End.XYComp > quant_E7.0[quant_breaks-1] & End.XYComp < quant_E7.0[quant_breaks] )
t.test(test_5$Avg.Density..cells.per.12.radii.,test_6$Avg.Density..cells.per.12.radii.)



