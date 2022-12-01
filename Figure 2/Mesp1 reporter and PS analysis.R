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
loadfonts()

library("RColorBrewer")
display.brewer.all()
brewer.pal(n = 9, name = "Paired")

#MaMuT_scale_factor = 3.942282
#microns_per_pixel = 0.380490284561

#==============
#labels and colors
#==============
x_labels <- c( "E6.5", "E6.75", "E7.0")
fills65 <-c("#64ABE3","#FDD8B5")
fills67 <-c("#15B2D1","#DFCE9D")
fills70 <-c("#65CBDA","#F9D199")
sand <- c(fills65[2],fills67[2],fills70[2])
water <- c(fills65[1],fills67[1],fills70[1])




#==============
# read data
#==============
E650_spots <-read_ods(path="E6.5 Mesp1lineage reporter initiation.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)
E675_spots <-read_ods(path="E6.75 Mesp1lineage reporter initiation.ods", sheet=1, formula_as_formula=FALSE, verbose=FALSE)
E650_landmarks <- read_ods(path="E6.5 primstreak and antmidline position.ods", sheet=1, row_names=TRUE, formula_as_formula=FALSE, verbose=FALSE)
E675_landmarks <- read_ods(path="E6.75 primstreak and antmidline position.ods", sheet=1, row_names=TRUE, formula_as_formula=FALSE, verbose=FALSE)


#==============
# assemble tracks
#==============
assemble_tracks <- function( spot_df ) { 
  track_list <- c()
  start_t <- c()
  end_t <- c()
  start_x <- c()
  end_x <- c()
  start_y <- c()
  end_y <- c()
  start_z <- c()
  end_z <- c()
  start_int <- c()
  end_int  <- c()
  
  for(i in seq(1,nrow(spot_df),2) ) {
    spot_1 <- spot_df[i,]
    spot_2 <- spot_df[i+1,]

    
    if ( spot_1$`Track ID` == spot_2$`Track ID` ) {
      #print ( paste0(spot_1$`Spot ID`,':',spot_2$`Spot ID` ) )
      if ( spot_1$T  > spot_2$T ) { #correct temporal order
        spot_3 = spot_2
        spot_2 = spot_1
        spot_1 = spot_3
      }

      track_list <- append(track_list, spot_1$`Track ID`)
      start_t <- append(start_t, spot_1$T)
      end_t <- append(end_t, spot_2$T)
      start_x <- append(start_x, spot_1$X)
      end_x <- append(end_x, spot_2$X)     
      start_y <- append(start_y, spot_1$Y)
      end_y <- append(end_y, spot_2$Y)
      start_z <- append(start_z, spot_1$Z)
      end_z <- append(end_z, spot_2$Z)  
      start_int <- append(start_int, spot_1$`Mean intensity ch1`)
      end_int <- append(end_int, spot_2$`Mean intensity ch1`)
      
      
    }

  }
  df <- data.frame( start_t, end_t, start_x, end_x, start_y, end_y, start_z, end_z, start_int, end_int, row.names=track_list)#, names=c("Begin T","End T","Begin X","End X","Begin Y","End Y","Begin Z","End Z","Begin Int","End Int") )
  names(df) <- c("Begin.T","End.T","Begin.X","End.X","Begin.Y","End.Y","Begin.Z","End.Z","Begin.Int","End.Int")
  return( df )
}
E650_tracks <- assemble_tracks(E650_spots)
E650_tracks$Stage <- "E6.5"
E675_tracks <- assemble_tracks(E675_spots)
E675_tracks$Stage <- "E6.75"



#==============
# calculate features
#==============
calculate_xz_features <- function( track_df, landmark_df ) { 
  as.character()
  
  
  track_df$Duration <- (track_df$End.T - track_df$Begin.T) * 6 #6 timeframes per minute
  
  track_df$DisplacementXZ <- sqrt((track_df$End.X-track_df$Begin.X)^2 + (track_df$End.Z-track_df$Begin.Z)^2)
  track_df$VelocityXZ <- track_df$DisplacementXZ * 60 / track_df$Duration
  
  track_df$Displacement <- sqrt((track_df$End.X-track_df$Begin.X)^2 + (track_df$End.Y-track_df$Begin.Y)^2 + (track_df$End.Z-track_df$Begin.Z)^2)
  track_df$Velocity <- track_df$Displacement * 60 / track_df$Duration  

  track_df$PS_Start.Distance <- sqrt((landmark_df[as.character(track_df$Begin.T), "Primitive Streak X"]-track_df$Begin.X)^2 + (landmark_df[as.character(track_df$Begin.T), "Primitive Streak Z"]-track_df$Begin.Z)^2)
  track_df$PS_End.Distance <- sqrt((landmark_df[as.character(track_df$End.T), "Primitive Streak X"]-track_df$End.X)^2 + (landmark_df[as.character(track_df$End.T), "Primitive Streak Z"]-track_df$End.Z)^2)
  track_df$AM_End.Distance <- sqrt((landmark_df[as.character(track_df$End.T), "Anterior Midline X"]-track_df$End.X)^2 + (landmark_df[as.character(track_df$End.T), "Anterior Midline Z"]-track_df$End.Z)^2)
  
  track_df$PS_to_Visible.EstTime <- track_df$PS_Start.Distance / (track_df$Velocity)
  track_df$PS_to_Tracked.EstTime <-  track_df$PS_to_Visible.EstTime + (track_df$Duration / 60)

  return( track_df )
}
E650_tracks <- calculate_xz_features(E650_tracks,E650_landmarks)
E675_tracks <- calculate_xz_features(E675_tracks,E675_landmarks)
All_tracks <- rbind(E650_tracks,E675_tracks)








#==============
#summary table plots
#==============
p1 <- All_tracks %>%
  ggplot( aes(x=Stage, y=Duration/60, fill=Stage, color=Stage, alpha=0.4)) +
  geom_violin( trim=FALSE, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels ) +#, limits = rev) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=water) +
  scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
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

p2 <- All_tracks %>%
  ggplot( aes(x=Stage, y=Velocity, fill=Stage, color=Stage, alpha=0.4)) +
  geom_violin( trim=FALSE, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(-10,250),breaks = seq(0, 200, len = 3)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels ) +#, limits = rev) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=water) +
  scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
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

p3 <- All_tracks %>%
  ggplot( aes(x=Stage, y=PS_to_Visible.EstTime, fill=Stage, color=Stage, alpha=0.4)) +
  geom_violin( trim=FALSE, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  scale_y_continuous(limits=c(-2,12),breaks = seq(0, 10, len = 3)) +
  scale_x_discrete(labels= x_labels ) +#, limits = rev) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=water) +
  scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
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

p4 <- All_tracks %>%
  ggplot( aes(x=Stage, y=PS_to_Tracked.EstTime, fill=Stage, color=Stage, alpha=0.4)) +
  geom_violin( trim=FALSE, lwd=2 ) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #coord_flip() + 
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete(labels= x_labels ) +#, limits = rev) +
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual(values=water) +
  scale_fill_manual(values=sand) +
  theme(
    legend.position="none",
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.y = element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
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

(p1 + theme(plot.margin = margin(0, 100, 0, 0)) ) | (p2 + theme(plot.margin = margin(0, 100, 0, 100)) ) | (p3 + theme(plot.margin = margin(0, 100, 0, 100)) ) | (p4 + theme(plot.margin = margin(0, 0, 0, 100)) )

# ( p1 / p2 ) | ( p3 / p4 ) 
# p1
# p2
# p3
# p4


#==============
#t-test
#==============
PosE65 <- All_tracks[which(All_tracks$Stage=="E6.5"),]
PosE67 <- All_tracks[which(All_tracks$Stage=="E6.75"),]

t.test(PosE65$Duration,PosE67$Duration)

t.test(PosE65$Velocity,PosE67$Velocity)

t.test(PosE65$PS_to_Visible.EstTime,PosE67$PS_to_Visible.EstTime)

t.test(PosE65$PS_to_Tracked.EstTime,PosE67$PS_to_Tracked.EstTime)
