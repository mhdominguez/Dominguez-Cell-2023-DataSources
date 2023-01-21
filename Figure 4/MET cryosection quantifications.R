#install.packages('readODS')
#font_import()
library(ggplot2)
library(readODS)
library(tidyverse)
library(viridis)
library(dplyr)
#library(extrafont)
library(patchwork)
#loadfonts()

library("RColorBrewer")
display.brewer.all()
brewer.pal(n = 9, name = "Paired")


#==============
#labels and colors
#==============
paired_colors <- rev(brewer.pal(name = "Dark2", n=6)[c(1,3)])



#==============
# read data and add factors
#==============
resultTable_cor <- read.csv(file = 'cryosection aCat-Pancad correlation.csv', sep=',')
resultTable_cor$Comparison <- as.factor(gsub("-Image[0-9]+","",resultTable_cor$Images))
resultTable_cor$Region <- gsub("[0-9]*","",resultTable_cor$Region)
resultTable_cor$Region <- sub("ect","Neuroepithelium",resultTable_cor$Region)
resultTable_cor$Region <- sub("end","Foregut Endoderm",resultTable_cor$Region)
resultTable_cor$Region <- sub("mes","MEF2C+ Mesoderm",resultTable_cor$Region)
resultTable_cor$Region <- as.factor(resultTable_cor$Region )
resultTable_cor$Stage <- sub("EHF","EHF",resultTable_cor$Stage)
resultTable_cor$Stage <- sub("Som","Early Somite",resultTable_cor$Stage)
resultTable_cor$Stage <- as.factor(resultTable_cor$Stage )

resultTable_intensity <- read.csv(file = 'cryosection aCat-Pancad intensity.csv', sep=',')
resultTable_intensity$Stage <- sub("EHF","EHF",resultTable_intensity$Stage)
resultTable_intensity$Stage <- sub("Som","Early Somite",resultTable_intensity$Stage)
resultTable_intensity$Stage <- as.factor(resultTable_intensity$Stage )
resultTable_intensity$Region <- gsub("Image[0-9]*:","",resultTable_intensity$Label)
resultTable_intensity$Region <- gsub("[^a-zA-Z]","",resultTable_intensity$Region)
resultTable_intensity$Region <- sub("ect","Neuroepithelium",resultTable_intensity$Region)
resultTable_intensity$Region <- sub("end","Foregut Endoderm",resultTable_intensity$Region)
resultTable_intensity$Region <- sub("mes","MEF2C+ Mesoderm",resultTable_intensity$Region)
resultTable_intensity$Region <- as.factor(resultTable_intensity$Region )




#==============
# Colocalization Pancad-aCatenin
#==============
Region.labels <- c("Foregut Endoderm"="Foregut\nEndoderm","MEF2C+ Mesoderm"="MEF2C+\nMesoderm","Neuroepithelium"="Neuro-\nepithelium")
p1 <- resultTable_cor %>% filter( Comparison == "C2 & C3" ) %>%
  ggplot( aes(x=Stage, y=Rtotal, color=Stage)) +
  #geom_violin( trim=FALSE, lwd=2 ) +
  #geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8, width=0.2) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25, show.legend=FALSE) +
  #coord_flip() + 
  facet_grid( ~Region,switch="both",labeller = labeller(Region = Region.labels))+
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  scale_y_continuous(trans="sqrt") +
  scale_x_discrete( limits = rev) + #labels= x_labels ) +#
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual (values=paired_colors) +
  #scale_fill_manual(values=paired_colors) +
  guides(color = guide_legend(reverse=TRUE)) +
  theme(
    #legend.position = c(0.16, 0.2),
    legend.position = c(0.5, 0.85),
    legend.box.background = element_rect(colour="#A0A0A0",size=1),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
    legend.title = element_text(size=56,vjust=1,hjust=0.5,color="black",family="Liberation Sans"),
    legend.text = element_text(size=48,color="black",family="Liberation Sans", angle=0, vjust=0.65, hjust=0),
    plot.title = element_text(size=56,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.y = element_text(size=48,color="black",family="Liberation Sans", angle=0, vjust=0.65, hjust=0.5),
    axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=48,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #panel.background = element_rect(fill="Light Gray"),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White"),
    strip.placement = "outside"
  ) +  
  xlab("")
ggsave( plot = p1, height=7, width=12, filename = "aCat-panCad localization.png", units = "in",dpi = 100 )

EHF <- resultTable_cor %>% filter( Comparison == "C2 & C3" & Stage == "EHF" )
Som <- resultTable_cor %>% filter( Comparison == "C2 & C3" & Stage == "Early Somite" )

t.test(EHF[which(EHF$Region=="Neuroepithelium"),]$Rtotal,Som[which(Som$Region=="Neuroepithelium"),]$Rtotal)
t.test(EHF[which(EHF$Region=="Foregut Endoderm"),]$Rtotal,Som[which(Som$Region=="Foregut Endoderm"),]$Rtotal)
t.test(EHF[which(EHF$Region=="MEF2C+ Mesoderm"),]$Rtotal,Som[which(Som$Region=="MEF2C+ Mesoderm"),]$Rtotal)

#t.test(EHF[which(EHF$Region=="Neuroepithelium"),]$Rcoloc,Som[which(Som$Region=="Neuroepithelium"),]$Rcoloc)
#t.test(EHF[which(EHF$Region=="Foregut Endoderm"),]$Rcoloc,Som[which(Som$Region=="Foregut Endoderm"),]$Rcoloc)
#t.test(EHF[which(EHF$Region=="MEF2C+ Mesoderm"),]$Rcoloc,Som[which(Som$Region=="MEF2C+ Mesoderm"),]$Rcoloc)



#==============
# Intensity of panCad, aCat, Mef2c
#==============
Ch.labels <- c("1"="MEF2C","2"="ɑ-CATENIN","3"="pan-\nCADHERIN")
p2 <- resultTable_intensity %>% filter( Region=="MEF2C+ Mesoderm" ) %>% filter( Ch == 1 | Ch == 2 | Ch == 3 ) %>% #Ch2 is aCat
  ggplot( aes(x=Stage, y=Median-Min, color=Stage)) +
  #geom_violin( trim=FALSE, lwd=2 ) +
  #geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8, width=0.18) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #coord_flip() + 
  facet_wrap( ~Ch,scales="free_y",switch="x",labeller = labeller(Ch = Ch.labels))+
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete( limits = rev) + #labels= x_labels ) +#
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual (values=paired_colors) +
  guides(color = guide_legend(reverse=TRUE)) +
  theme(
    legend.position = "none",
    legend.box.background = element_rect(colour="#A0A0A0",size=1),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
    legend.title = element_text(size=56,vjust=1,hjust=0.5,color="black",family="Liberation Sans"),
    legend.text = element_text(size=48,color="black",family="Liberation Sans", angle=0, vjust=0.65, hjust=0),
    plot.title = element_text(size=56,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.y = element_blank(), #element_text(size=48,color="black",family="Liberation Sans", angle=0, vjust=0.65, hjust=0.5),
    axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=48,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(3.5, "lines"),
    #panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White"),
    strip.placement = "outside"
  ) +  
  xlab("")
ggsave( plot = p2, height=7, width=12, filename = "Mef2c-aCat-panCad-meso intensity.png", units = "in",dpi = 100 )

EHF <- resultTable_intensity %>% filter( Stage == "EHF" & Region=="MEF2C+ Mesoderm")
Som <- resultTable_intensity %>% filter( Stage == "Early Somite" & Region=="MEF2C+ Mesoderm")

for ( i in seq(1,3) ) {
  print ( t.test(EHF[which(EHF$Ch==i),]$Median-EHF[which(EHF$Ch==i),]$Min,Som[which(Som$Ch==i),]$Median-Som[which(Som$Ch==i),]$Min) )
}









#==============
# read data and add factors
#==============
resultTable_cor <- read.csv(file = 'cryosection Jup-nCad correlation.csv', sep=',')
resultTable_cor$Comparison <- as.factor(gsub("-Image[0-9]+","",resultTable_cor$Images))
resultTable_cor$Region <- gsub("[0-9]*","",resultTable_cor$Region)
resultTable_cor$Region <- sub("ect","Neuroepithelium",resultTable_cor$Region)
resultTable_cor$Region <- sub("end","Foregut Endoderm",resultTable_cor$Region)
resultTable_cor$Region <- sub("mes","MEF2C+ Mesoderm",resultTable_cor$Region)
resultTable_cor$Region <- as.factor(resultTable_cor$Region )
resultTable_cor$Stage <- sub("EHF","EHF",resultTable_cor$Stage)
resultTable_cor$Stage <- sub("Som","Early Somite",resultTable_cor$Stage)
resultTable_cor$Stage <- as.factor(resultTable_cor$Stage )

resultTable_intensity <- read.csv(file = 'cryosection Jup-nCad intensity.csv', sep=',')
resultTable_intensity$Stage <- sub("EHF","EHF",resultTable_intensity$Stage)
resultTable_intensity$Stage <- sub("Som","Early Somite",resultTable_intensity$Stage)
resultTable_intensity$Stage <- as.factor(resultTable_intensity$Stage )
resultTable_intensity$Region <- gsub("Image[0-9]*:","",resultTable_intensity$Label)
resultTable_intensity$Region <- gsub("[^a-zA-Z]","",resultTable_intensity$Region)
resultTable_intensity$Region <- sub("ect","Neuroepithelium",resultTable_intensity$Region)
resultTable_intensity$Region <- sub("end","Foregut Endoderm",resultTable_intensity$Region)
resultTable_intensity$Region <- sub("mes","MEF2C+ Mesoderm",resultTable_intensity$Region)
resultTable_intensity$Region <- as.factor(resultTable_intensity$Region )













#==============
# Colocalization gCat-nCad
#==============
Region.labels <- c("Foregut Endoderm"="Foregut\nEndoderm","MEF2C+ Mesoderm"="MEF2C+\nMesoderm","Neuroepithelium"="Neuro-\nepithelium")
p1 <- resultTable_cor %>% filter( Region == "MEF2C+ Mesoderm" & Comparison == "C2 & C3" ) %>%
  ggplot( aes(x=Stage, y=Rtotal, color=Stage)) +
  #geom_violin( trim=FALSE, lwd=2 ) +
  #geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8, width=0.2) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25, show.legend=FALSE) +
  #coord_flip() + 
  facet_grid( ~Region,switch="both",labeller = labeller(Region = Region.labels))+
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete( limits = rev) + #labels= x_labels ) +#
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual (values=paired_colors) +
  #scale_fill_manual(values=paired_colors) +
  guides(color = guide_legend(reverse=TRUE)) +
  theme(
    legend.position = "none",
    #legend.position = c(0.8, 0.4),
    legend.box.background = element_rect(colour="#A0A0A0",size=1),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
    legend.title = element_text(size=56,vjust=1,hjust=0.5,color="black",family="Liberation Sans"),
    legend.text = element_text(size=48,color="black",family="Liberation Sans", angle=0, vjust=0.65, hjust=0),
    plot.title = element_text(size=56,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.y = element_text(size=48,color="black",family="Liberation Sans", angle=0, vjust=0.65, hjust=0.5),
    axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=48,color="black",family="Liberation Sans"),
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
ggsave( plot = p1, height=7, width=4.7, filename = "gCat-nCad localization.png", units = "in",dpi = 100 )

EHF <- resultTable_cor %>% filter( Comparison == "C2 & C3" & Stage == "EHF" )
Som <- resultTable_cor %>% filter( Comparison == "C2 & C3" & Stage == "Early Somite" )

t.test(EHF[which(EHF$Region=="MEF2C+ Mesoderm"),]$Rtotal,Som[which(Som$Region=="MEF2C+ Mesoderm"),]$Rtotal)
#t.test(EHF[which(EHF$Region=="MEF2C+ Mesoderm"),]$Rcoloc,Som[which(Som$Region=="MEF2C+ Mesoderm"),]$Rcoloc)



#==============
# Intensity of nCad, gCat
#==============
Ch.labels <- c("1"="DAPI","2"="N-CADHERIN","3"="γ-CATENIN","4"="MEF2C")
p2 <- resultTable_intensity %>% filter( Region=="MEF2C+ Mesoderm" ) %>% filter( Ch == 2 | Ch == 3 ) %>%
  ggplot( aes(x=Stage, y=Median-Min, color=Stage)) +
  #geom_violin( trim=FALSE, lwd=2 ) +
  #geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(color=Stage), size=8, alpha=0.8, width=0.18) +
  geom_crossbar(stat = "summary", fun=mean,color="black",size=1.25) +
  #coord_flip() + 
  facet_wrap( ~Ch,scales="free_y",switch="x",labeller = labeller(Ch = Ch.labels))+
  #facet_wrap( ncol = 1, vars(Axis), labeller = labeller(Axis = new_xy_labels) )+
  #ylab("fraction of apical-basal length") +
  #ggtitle(paste(expression(italic("Rosa26"^{"Ai66"})), " per Heart"))) +
  #ggtitle("basic linear ratios, by heart") + 
  theme_light() +
  #scale_y_continuous(limits=c(0.25,0.75)) +
  #scale_y_continuous(trans="sqrt") +
  scale_x_discrete( limits = rev) + #labels= x_labels ) +#
  #scale_color_brewer(palette="YlOrRd") +
  #scale_fill_brewer(palette="YlOrRd") +
  #scale_color_manual(values=c("black", "black", "black")) +
  scale_color_manual (values=paired_colors) +
  guides(color = guide_legend(reverse=TRUE)) +
  theme(
    legend.position = "none",
    legend.box.background = element_rect(colour="#A0A0A0",size=1),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
    legend.title = element_text(size=56,vjust=1,hjust=0.5,color="black",family="Liberation Sans"),
    legend.text = element_text(size=48,color="black",family="Liberation Sans", angle=0, vjust=0.65, hjust=0),
    plot.title = element_text(size=56,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title = element_blank(),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    axis.ticks = element_blank(), #element_line(color="#444444",size=3),
    #axis.text.x = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    axis.text.y = element_blank(), #element_text(size=48,color="black",family="Liberation Sans", angle=0, vjust=0.65, hjust=0.5),
    axis.text.x = element_blank(), #element_text(size=60,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=48,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(4, "lines"),
    #panel.grid.major.y = element_line(colour="#A0A0A0",size=1),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(10, 30, 0, 0),
    rect = element_blank(),
    #panel.background = element_rect(fill=fills65[2]),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill="White"),
    strip.placement = "outside"
  ) +  
  xlab("")
ggsave( plot = p2, height=7, width=8.7, filename = "nCad-gCat-meso intensity.png", units = "in",dpi = 100 )

EHF <- resultTable_intensity %>% filter( Stage == "EHF" & Region=="MEF2C+ Mesoderm")
Som <- resultTable_intensity %>% filter( Stage == "Early Somite" & Region=="MEF2C+ Mesoderm")

for ( i in seq(2,3) ) {
  print ( t.test(EHF[which(EHF$Ch==i),]$Median-EHF[which(EHF$Ch==i),]$Min,Som[which(Som$Ch==i),]$Median-Som[which(Som$Ch==i),]$Min) )
}





