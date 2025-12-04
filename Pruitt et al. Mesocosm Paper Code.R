rm(list=ls())
setwd("")
getwd()
#Packages needed 
library(chron) #chron helps deal with dates and times in R
library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggthemes)
library(ggfortify)
library(gridExtra)
library(RColorBrewer)
library(pracma)
library(patchwork)
library(FSA)
library(car)
library(fmsb)
library(ggridges)
library(viridis)
library(plyr)
library(rlang)
library(lme4)
library(ggResidpanel)
library(emmeans)
library(rstatix)
library(ggpubr)
library(influence.ME)
library(ggh4x)
library(corrplot)
library(lme4)
library(nlme)
library(TMB)
library(glmmTMB)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

#dev.off()
#set strings as factors to false
options(stringsAsFactors = FALSE)

#########Chla#######
data<-read.csv(file="mesocc.chla.5.7.24.csv")
#Average for total
dt <- data.table(data)

z_scores <- scale(abs(dt$Total))
outliers<-which(abs(z_scores) > 3)
dtnew <- dt[-outliers, ]

# Summarize data including Cyano, Green, and Diatom
dt.out <- dtnew[, .(
  mean_Total = mean(Total, na.rm = TRUE),
  sd_Total = sd(Total, na.rm = TRUE),
  count_Total = .N,
  mean_Cyano = mean(Cyano, na.rm = TRUE),
  sd_Cyano = sd(Cyano, na.rm = TRUE),
  mean_Green = mean(Green, na.rm = TRUE),
  sd_Green = sd(Green, na.rm = TRUE),
  mean_Diatom = mean(Diatom, na.rm = TRUE),
  sd_Diatom = sd(Diatom, na.rm = TRUE)
), by = c("substrate", "temp", "substrate.temp", "day","Day")]

# Calculate standard errors
dt.out[, `:=`(
  se_Total = sd_Total / sqrt(count_Total),
  se_Cyano = sd_Cyano / sqrt(count_Total),
  se_Green = sd_Green / sqrt(count_Total),
  se_Diatom = sd_Diatom / sqrt(count_Total)
)]

# Reorder factors
dt.out$substrate.temp <- factor(dt.out$substrate.temp, levels = c("Cobble-20", "Leaf-20", "Cobble-23", "Leaf-23", "Cobble-26", "Leaf-26"))
dt.out$temp <- factor(dt.out$temp, levels = c("20°","23°", "26°"))
dt.out$Day <- factor(dt.out$Day, levels = c("Day 1", "Day 7", "Day 14", "Day 21"))
cobchla<-subset(dt.out, substrate=="Cobble")

chlafig<-ggplot(cobchla, aes(x = day, y =mean_Total,  group = temp,fill = temp)) +
  geom_line(
    aes(group = interaction(temp)),  # Group by temp and subs
    color = "black",
    size = 0.5 ) +
  geom_point(size = 4, color = "black",shape=22) +
  geom_errorbar(aes(ymax = (mean_Total+se_Total), ymin = (mean_Total-se_Total)), 
                stat="identity", width = 0.2)+
  scale_fill_manual(values = c("yellow", "orange", "red")) +
  scale_shape_manual(values = c(22, 23)) +#21, 
  xlab(expression("Day")) +
  ylab(expression(Total~Chla~(ug~cm^{"-2"}))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(
    axis.text.x = element_text(size = 12, color = "black", vjust = 1.5),
    axis.title.y = element_text(size = 13, color = "black", vjust = 1.5),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text.y=element_blank(),
    strip.background.y=element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title.x = element_text(size = 13, color = "black", vjust = 1.5),
    legend.position = "none"
  )+
  scale_y_continuous(limits=c(0,1.5))+
  scale_x_continuous(breaks=c(1,7,14,21),labels = c(1,7,14,21))
chlafig

df_longcob <- cobchla %>%
  pivot_longer(
    cols = c(mean_Cyano, mean_Green, mean_Diatom),
    names_to = "Chla_Type",
    values_to = "Chla_Value"
  )

SFig3<-ggplot(df_longcob, aes(x = factor(temp), y = Chla_Value, fill = Chla_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(. ~ Day) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0),limits=c(0,1.5))+
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+scale_fill_manual(values = c("cyan4", "burlywood4", "darkgreen")) +  # Custom fill colors
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.title.y = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12)
  ) +
  xlab(expression("Temperature (°C)")) +
  ylab(expression(Chla~(ug~cm^{"-2"}))) 

SFig3

#####CR####
data<-read.csv(file="mesocc.metab.forr.csv")
filtered_data <- data[!is.na(data$ER.mgm2h), ]
data<-filtered_data

#Average for CR
z_scores <- scale(abs(data$ER.mgm2h))
which(abs(z_scores) > 3)  #NONE
outliers<-which(abs(z_scores) > 3)

dt<-data.table(data)
dt.out<-dt[,list(mean=mean(ER.mgm2h, na.rm=TRUE),sd=sd(ER.mgm2h, na.rm=TRUE), count=length(ER.mgm2h)),
           by=c("sub.temp","Day","subs","temp","day")]
dt.out$se<-dt.out$sd/sqrt(dt.out$count)
dt.out$substrate.temp<-factor(dt.out$sub.temp, levels = c("Cobble-20", "Leaf-20", "Cobble-23","Leaf-23","Cobble-26","Leaf-26"))
dt.out$Day<-factor(dt.out$Day, levels = c("Day 1","Day 7","Day 14","Day 21"))
dt.out$temp<-factor(dt.out$temp, levels = c("20°","23°","26°"))

leafcr<-subset(dt.out, subs=="Leaves")
crplotall<-ggplot(leafcr, aes(x = day, y = abs(mean),  group = temp,fill = temp)) +
  geom_line(
    color = "black",
    size = 0.5 ) +
  geom_point(size = 4, color = "black",shape=23) +
  geom_errorbar(aes(ymax = abs(mean+se), ymin = abs(mean-se)), 
                stat="identity", width = 0.2)+
  scale_fill_manual(values = c("yellow", "orange", "red")) +
  xlab(expression("Day")) +
  ylab(expression("|CR|"~~(mg~O["2"]~m^{-2}~hr^{-1}))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(
    axis.title.x = element_text(size = 13, color = "black", vjust = 1.5),
    axis.title.y = element_text(size = 13, color = "black", vjust = 1.5),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text.y=element_blank(),
    strip.background.y=element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    legend.position = "none"
  )+
  scale_x_continuous(breaks=c(1,7,14,21),labels = c(1,7,14,21))
crplotall

####Penetrometer####
data<-read.csv(file="mesocc.penetrometer.12.6.23.csv",fileEncoding = "UTF-8")
dt<-data.table(data)
dt.out<-dt[,list(mean=mean(pen.kg, na.rm=TRUE),sd=sd(pen.kg, na.rm=TRUE), count=length(pen.kg)),
           by=c("temp","Day","day")]
dt.out$se<-dt.out$sd/sqrt(dt.out$count)
dt.out$temp<-factor(dt.out$temp, levels = c("20°", "23°", "26°"))

penfig<-ggplot(dt.out, aes(x = day, y =mean,  group = temp,fill = temp)) +
  geom_line(
    aes(group = interaction(temp)),  # Group by temp and subs
    color = "black",
    size = 0.5 ) +
  geom_point(size = 4, color = "black",shape=23) +
  geom_errorbar(aes(ymax = (mean+se), ymin = (mean-se)), 
                stat="identity", width = 0.2)+
  scale_fill_manual(values = c("yellow", "orange", "red")) +
  xlab(expression("Day")) +
  ylab(expression(Penetrometer~(kg))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 13, color = "black", vjust = 1.5),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text.y=element_blank(),
    strip.background.y=element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title.x = element_text(size = 13, color = "black", vjust = 1.5),
    legend.position = "none"
  )+
  scale_y_continuous(limits=c(.5,1.5))+
  scale_x_continuous(breaks=c(1,7,14,21),labels = c(1,7,14,21))
penfig

#Figure 4
ancillaryplot<-crplotall/penfig/chlafig
ancillaryplot

#####K and U####
data<-read.csv(file="K_modeloutputsandcalcsforpaper.csv",fileEncoding = "latin1")
data$temp  = factor(data$temp, levels=c("20°", "23°", "26°"))
data$day  = factor(data$day, levels=c("Day 7", "Day 14", "Day 21"))
cobble<-subset(data,substrate=="cobble")
leaf<-subset(data,substrate=="leaf")

#Fig3
cobblek<-ggplot() +
  geom_line(
    data = cobble,
    aes(
      y = k.hr,
      x = temp,
      group = interaction(nutrient, day)),color = "black",
    size = 0.5) +
  geom_pointrange(
    data = cobble,
    colour = "black",
    size = 1,  # Line thickness
    point.size = 4,  # Point size
    width = 0.8,
    shape = 21,
    aes(y = k.hr, x = temp, ymin = abs(k.low.hr), ymax = abs(k.high.hr), fill = temp)
  ) +
  scale_fill_manual(values = c("yellow", "orange", "red")) +  # Custom fill colors
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.title.y = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text.y=element_blank(),
    strip.background.y=element_blank())+
  xlab(expression("Temperature Treatment (°C)")) +
  ylab(expression(k~(hr^{"-1"}))) +
  scale_y_continuous(labels = scales::comma)+
  facet_grid(nutrient ~ day, scales = "free_y", axis.labels = "margins", axes = "all")

leafk<-ggplot() +
  geom_line(
    data = leaf,
    aes(
      y = k.hr,
      x = temp,
      group = interaction(nutrient, day)),color = "black",
    size = 0.5) +
  geom_pointrange(
    data = leaf,
    colour = "black",
    size = 1,  # Line thickness
    point.size = 4,  # Point size
    width = 0.8,
    shape = 21,
    aes(y = k.hr, x = temp, ymin = abs(k.low.hr), ymax = abs(k.high.hr), fill = temp)
  ) +
  scale_fill_manual(values = c("yellow", "orange", "red")) +  # Custom fill colors
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  xlab(expression("Temperature Treatment (°C)")) +
  ylab(expression(k~(hr^{"-1"}))) +
  #scale_y_continuous(labels = scales::comma,breaks=seq(0,max(leaf$k.high.hr),length.out=5))+
  facet_grid(nutrient ~ day, scales = "free_y", axis.labels = "margins", axes = "all")

cobblek|leafk

#SFig2
cobbleu<-ggplot() +
  geom_line(
    data = cobble,
    aes(
      y = U,
      x = temp,
      group = interaction(nutrient, day)),color = "black",
    size = 0.5) +
  geom_pointrange(
    data = cobble,
    colour = "black",
    size = 1,  # Line thickness
    point.size = 3,  # Point size
    width = 0.8,
    aes(y = U, x = temp, ymin = abs(U.low), ymax = abs(U.high), fill = temp, shape = day)
  ) +
  scale_fill_manual(values = c("yellow", "orange", "red")) +  # Custom fill colors
  theme_classic() +
  
  scale_shape_manual(values=c(21,22,23))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.title.y = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text.y=element_blank(),
    strip.background.y=element_blank())+
  xlab(expression("Temperature Treatment (°C)")) +
  ylab(expression(U~~(mg~m^{"-2"}~d^{"-1"})))+
  scale_y_continuous(labels = scales::comma)+
  facet_grid(nutrient ~ ., scales = "free_y", axis.labels = "margins", axes = "all")

leafu<-ggplot() +
  geom_line(
    data = leaf,
    aes(
      y = U,
      x = temp,
      group = interaction(nutrient, day)),color = "black",
    size = 0.5) +
  geom_pointrange(
    data = leaf,
    colour = "black",
    size = 1,
    aes(y = U, x = temp, ymin = abs(U.low), ymax = abs(U.high), fill = temp,shape=day)
  ) +
  scale_fill_manual(values = c("yellow", "orange", "red")) +  # Custom fill colors
  theme_classic() +
  scale_shape_manual(values=c(21,22,23))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.title.y = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text.y=element_blank(),
    strip.background.y=element_blank())+
  xlab(expression("Temperature Treatment (°C)")) +
  ylab(expression(U~~(mg~m^{"-2"}~d^{"-1"})))+
  facet_grid(nutrient ~ ., scales = "free_y", axis.labels = "margins", axes = "all")

cobbleuchla<-ggplot() +
  geom_line(
    data = cobble,
    aes(
      y = U.chla.mean.mgmgchlad,
      x = temp,
      group = interaction(nutrient, day)),color = "black",
    size = 0.5) +
  geom_pointrange(
    data = cobble,
    colour = "black",
    size = 1,  # Line thickness
    point.size = 4,  # Point size
    width = 0.8,
    aes(y = U.chla.mean.mgmgchlad, x = temp, ymin = abs(U.chla.low.mgmgchlad), ymax = abs(U.chla.high.mgmgchlad), fill = temp, shape = day)
  ) +
  scale_fill_manual(values = c("yellow", "orange", "red")) +  # Custom fill colors
  theme_classic() +
  scale_shape_manual(values=c(21,22,23))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.title.y = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text.y=element_blank(),
    strip.background.y=element_blank())+
  xlab(expression("Temperature Treatment (°C)")) +
  ylab(expression(U~Efficiency~(mg~"\u00B7"~mg~chla^{"-1"}~"\u00B7"~d^{"-1"})))+
  ylab(expression(U~Efficiency~(mg~mg~chla^{"-1"}~d^{"-1"})))+
  scale_y_continuous(labels = scales::comma)+
  facet_grid(nutrient ~ ., scales = "free_y", axis.labels = "margins", axes = "all")

leafudm<-ggplot() +
  geom_line(
    data = leaf,
    aes(
      y = U.dm.mean.mggdmd,
      x = temp,
      group = interaction(nutrient, day)),color = "black",
    size = 0.5) +
  geom_pointrange(
    data = leaf,
    colour = "black",
    size = 1,
    aes(y = U.dm.mean.mggdmd, x = temp, ymin = abs(U.dm.low.mggdmd), ymax = abs(U.dm.high.mggdmd), fill = temp,shape = day)
  ) +
  scale_fill_manual(values = c("yellow", "orange", "red")) +  # Custom fill colors
  theme_classic() +
  scale_y_continuous(labels = scales::label_number()) +
  scale_shape_manual(values=c(21,22,23))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.title.y = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  xlab(expression("Temperature Treatment (°C)")) +
ylab(expression(U~Efficiency~(mg~"\u00B7"~g~DM^{"-1"}~"\u00B7"~d^{"-1"})))+
  facet_grid(nutrient ~ ., scales = "free_y", axis.labels = "margins", axes = "all")


cobbleu|cobbleuchla|leafu|leafudm

data<-read.csv(file="K_modeloutputsandcalcs_simp.csv",fileEncoding = "latin1")
data$temp  = factor(data$temp, levels=c("20°", "23°", "26°"))
data$day  = factor(data$day, levels=c("Day 7", "Day 14", "Day 21"))

#Fig 6
nh4k<-ggplot(data, aes(x = bkd.nh4.ugl, y = nh4.k, group = temp, fill = temp, shape = day)) +
  geom_point(size = 5, color = "black") +
  # geom_smooth(aes(group = s), method = "lm", formula = y ~ x, se = FALSE,
  #             color = "black", linetype = "solid") +
  # geom_smooth(data = subset(data, s == "Cobble"),
  #             aes(group=s),  # no group/fill/shape here!
  #             method = "lm", formula = y ~ x, se = FALSE,
  #             color = "black", linetype = "solid") +
  scale_fill_manual(values = c("yellow", "orange", "red")) +
  scale_shape_manual(values = c(21, 22, 24)) +  # Customize shapes if needed
  xlab(expression("Bkd NH"[4]*phantom()^"+"*"-N ("*mu*"g L"^-1*")"))+
  ylab(expression("Mean NH"[4]*phantom()^"+"*"-N k (hr"^{-1}*")"))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 13, color = "black", vjust = 1.5),
    axis.title.x = element_text(size = 13, color = "black", vjust = 1.5),
    strip.text.y = element_blank(),
    strip.background.y = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    legend.position = "none"
  ) +
  facet_grid(. ~ s, scales = "free_x")
srpk<-ggplot(data, aes(x = bkd.srp.ugl, y = srp.k, group = temp, fill = temp, shape = day)) +
  geom_point(size = 5, color = "black") +
  # geom_smooth(aes(group = s), method = "lm", formula = y ~ x, se = FALSE,
  #             color = "black", linetype = "solid") +
  geom_smooth(data = data,
              aes(group=s),  # no group/fill/shape here!
              method = "lm", formula = y ~ x, se = FALSE,
              color = "black", linetype = "solid") +
  scale_fill_manual(values = c("yellow", "orange", "red")) +
  scale_shape_manual(values = c(21, 22, 24)) +  # Customize shapes if needed
  xlab(expression("Bkd SRP ("*mu*"g L"^-1*")"))+
  ylab(expression("Mean SRP "*k*" (hr"^{-1}*")"))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 13, color = "black", vjust = 1.5),
    axis.title.x = element_text(size = 13, color = "black", vjust = 1.5),
    strip.text.x = element_blank(),
    strip.background.x = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    legend.position = "none"
  ) +
  facet_grid(. ~ s, scales = "free_x")
nh4k/srpk

SFig4<-ggplot(data, aes(x = nh4.k, y = srp.k, group = temp, fill = temp, shape = day)) +
  geom_point(size = 5, color = "black") +
  # geom_smooth(aes(group = s), method = "lm", formula = y ~ x, se = FALSE,
  #             color = "black", linetype = "solid") +
  geom_smooth(data = subset(data, s == "Cobble"),
              aes(group=s),  # no group/fill/shape here!
              method = "lm", formula = y ~ x, se = FALSE,
              color = "black", linetype = "solid") +
  scale_fill_manual(values = c("yellow", "orange", "red")) +
  scale_shape_manual(values = c(21, 22, 24)) +  # Customize shapes if needed
  xlab(expression("Mean NH"[4]*phantom()^"+"*"-N k (hr"^{-1}*")"))+
  ylab(expression("Mean SRP "*k*" (hr"^{-1}*")"))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 13, color = "black", vjust = 1.5),
    axis.title.x = element_text(size = 13, color = "black", vjust = 1.5),
    strip.text.x = element_blank(),
    strip.background.x = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    legend.position = "none"
  ) +
  theme(plot.margin = margin(.5,.5,.5,.5, "cm"))+
  facet_grid(. ~ s, scales = "free_x")
SFig4

######Figure 2######
Mesos<-read.csv("mesoccuptakeregrs.csv")

M_All.mean<-data.table(Mesos)
#Reset orders
M_All.mean$day = factor(M_All.mean$day, levels=c("Day 7", "Day 14","Day 21"))
M_All.mean$temperature = factor(M_All.mean$temperature, levels=c("20°", "23°","26°"))
M_All.mean$substrate = factor(M_All.mean$substrate, levels=c("Leaf", "Cobble"))
#Scale by first time point
filtered_data <- M_All.mean %>%
  filter(time.n == 1)
filtered_data<-filtered_data[,list(meanfirst=mean(ln.t, na.rm=TRUE)),
                             by=c("mesocosm", "nut","day")]
M_All.mean <- M_All.mean %>%
  left_join(filtered_data, by = c("mesocosm", "nut","day"))

M_All.mean <- M_All.mean %>%
  mutate(
    mean_divided = ln.t / meanfirst  # Divide 'mean' by the previously calculated 'mean_value'
  )

#subset
nh4<-subset(M_All.mean,nut=="NH4-N")
leafnh4<-subset(nh4,substrate=="Leaf")
no3<-subset(M_All.mean,nut=="NO3-N")
leafno3<-subset(no3,substrate=="Leaf")
srp<-subset(M_All.mean,nut=="SRP")
leafsrp<-subset(srp,substrate=="Leaf")

cobnh4<-subset(nh4,substrate=="Cobble")
cobno3<-subset(no3,substrate=="Cobble")
cobsrp<-subset(srp,substrate=="Cobble")

#plot it

#Updated code
nh4.meso.cob <- ggplot(cobnh4, aes(x = time, y = mean_divided, fill = as.factor(temperature), group = meso.number)) +
  geom_smooth(data = subset(cobnh4, 
                            (meso.number == "1" | meso.number == "2" | meso.number == "3" | meso.number == "4") & sig == "y"),
              aes(group = meso.number), method = 'lm', se = FALSE, color = "black", linetype = "solid", size = 0.7) +  # Black and thinner line
  #ylab(expression(Ln(Tracer~NH4~N))) +
  xlab(expression(Time~(hr))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  geom_point(size = 3, aes(shape = factor(meso.number), fill = factor(temperature)), color = "black", stroke = 0.5) +
  scale_color_manual(values = c("yellow", "orange", "red")) +
  scale_fill_manual(values = c("yellow", "orange", "red")) +
  scale_y_continuous(limits = c(-.7, 1.1)) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +  # Assign unique shapes for each mesocosm
  theme(strip.text.x = element_text(size = 13, color = "black"),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.background.y = element_blank()) + 
  
  facet_grid(day~temperature)

nh4.meso.cob

srp.meso.cob <- ggplot(cobsrp, aes(x = time, y = mean_divided, fill = as.factor(temperature), group = meso.number)) +
  geom_smooth(data = subset(cobsrp, 
                            (meso.number == "1" | meso.number == "2" | meso.number == "3" | meso.number == "4") & sig == "y"),
              aes(group = meso.number), method = 'lm', se = FALSE, color = "black", linetype = "solid", size = 0.7) +  # Black and thinner line
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(Ln(Tracer~SRP))) +
  xlab(expression(Time~(hr))) +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") + #geom_point(size = 3, aes(shape = factor(meso.number), color = factor(temperature))) +
  geom_point(size = 3, aes(shape = factor(meso.number), fill = factor(temperature)), color = "black", stroke = 0.5) +
  scale_color_manual(values = c("yellow", "orange", "red")) +
  scale_fill_manual(values = c("yellow", "orange", "red")) +
  scale_y_continuous(limits = c(-.7, 1.1)) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +# Assign unique shapes for each mesocosm
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_blank()) + 
  facet_grid(day~temperature)

srp.meso.cob

cobblemesosall<-nh4.meso.cob/srp.meso.cob

nh4.meso.leaf <- ggplot(leafnh4, aes(x = time, y = mean_divided, fill = as.factor(temperature), group = meso.number)) +
  geom_smooth(data = subset(leafnh4, 
                            (meso.number == "1" | meso.number == "2" | meso.number == "3" | meso.number == "4") & sig == "y"),
              aes(group = meso.number), method = 'lm', se = FALSE, color = "black", linetype = "solid", size = 0.7) +  # Black and thinner line
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(Ln(Tracer~NH4~N))) +
  xlab(expression(Time~(hr))) +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  geom_point(size = 3, aes(shape = factor(meso.number), fill = factor(temperature)), color = "black", stroke = 0.5) +
  scale_color_manual(values = c("yellow", "orange", "red")) +
  scale_fill_manual(values = c("yellow", "orange", "red")) +
  scale_y_continuous(limits = c(-2, 1.1)) +
  scale_shape_manual(values = c(21, 22, 23, 24)) + # Assign unique shapes for each mesocosm
  theme(strip.text.x = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.y = element_text(size = 13, color = "black"))+
  facet_grid(day~temperature)
nh4.meso.leaf

srp.meso.leaf <- ggplot(leafsrp, aes(x = time, y = mean_divided, fill = as.factor(temperature), group = meso.number)) +
  geom_smooth(data = subset(leafsrp, 
                            (meso.number == "1" | meso.number == "2" | meso.number == "3" | meso.number == "4") & sig == "y"),
              aes(group = meso.number), method = 'lm', se = FALSE, color = "black", linetype = "solid", size = 0.7) +  # Black and thinner line
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(Ln(Tracer~SRP))) +
  xlab(expression(Time~(hr))) +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  geom_point(size = 3, aes(shape = factor(meso.number), fill = factor(temperature)), color = "black", stroke = 0.5) +
  scale_color_manual(values = c("yellow", "orange", "red")) +
  scale_fill_manual(values = c("yellow", "orange", "red")) +
  scale_y_continuous(limits = c(-.1, 1.1)) +
  scale_shape_manual(values = c(21, 22, 23, 24)) + # Assign unique shapes for each mesocosm
  theme(strip.text.x = element_blank(),
        strip.text.y = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.background.x=element_blank())+
  facet_grid(day~temperature)
srp.meso.leaf

leafmesosall<-nh4.meso.leaf/srp.meso.leaf
cobblemesosall|leafmesosall

######PCA LEAF####
rm(list=ls())
setwd("")
getwd()

list.of.packages<-c("ggplot2", "dplyr", "gtools", "stats", "lme4", "car", "arm",
                    "multcomp","lattice","performance", "MASS", "buildmer", "mgcv", "caret", 
                    "tidyverse", "vegan", "VIM", "reshape2", "ggrepel", "Hmisc")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)>0) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

# Load data
df <- read.csv("mesocc.pca.leaf.csv")

df$temp <- factor(df$temp, levels=c("20", "23", "26"))
df$subs  = factor(df$subs, levels=c("Leaf"))

#df$season <- as.character(df$season)

# Subset environmental variables (for PCA)
df_env  = subset(df , select = c(ER.mgm2h,chla.ugcm2,nh4.bkd,srp.bkd,no3.bkd,doc.mgl,pen.g))

# Run scaled PCA
df_PCA <- rda(df_env, scale = TRUE)
summary(eigenvals(df_PCA))

# Bind PCs to data
df.corr <- cbind.data.frame(df,
                            scores(df_PCA, display = "sites", scaling = 1))


# pull out correlation variables
data.correlation = subset(df.corr, select = c(Day,Date,Substrate,sub.day,sub.temp,subs,temp,bank,nh4.k,srp.k,
                                              ER.mgm2h,chla.ugcm2,nh4.bkd,srp.bkd,no3.bkd,doc.mgl,pen.g,
                                              PC1,PC2))

correlation <- rcorr(as.matrix(data.correlation[,9:19]), type="pearson")

correlation


#Repeat as necessary
# Extract response variable from dataframe
dat <- cbind.data.frame(response = df$nh4.k,
                        scores(df_PCA, display = "sites", scaling = 1))

# Fit response to PCA
mod <- lm(response ~ PC1 + PC2, data = dat)
summary(mod) #PC1 and 2 signif
mod <- lm(response ~ PC2, data = dat)
summary(mod)

# Extract env and response scores and make tables (this don't need any changes)
Resp_sco <- data.frame(summary(df_PCA)$species[,1:3]) # get the species PC1, PC2 & PC3 scores
Site_sco <- data.frame(summary(df_PCA)$sites[,1:3])
Env_sco <- data.frame(summary(df_PCA)$biplot[,1:3])
Resp_tbl <- as_tibble(Resp_sco)
Site_tbl <- as_tibble(Site_sco)
Env_tbl <- as_tibble(Env_sco)
Resp_tbl <- mutate(Resp_tbl, vgntxt=rownames(Resp_sco),
                   ccatype = "species")
Site_tbl <- mutate(Site_tbl, vgntxt=rownames(Site_tbl),
                   ccatype = "sites")
Env_tbl <- mutate(Env_tbl, vgntxt=rownames(Env_sco),
                  ccatype = "bp")

# Use same env variables as in PCA
Env_scores <- Resp_tbl%>% filter(vgntxt %in% c("ER.mgm2h","chla.ugcm2","nh4.bkd","srp.bkd","no3.bkd","doc.mgl","pen.g"))

# Use same response variable
Response_scores <- Resp_tbl%>% filter(vgntxt %in% c("srp.k", "nh4.k")) # Same as response variable

# Extract grouping factor i.e., day/night. This is used for coloring points
df_group<-read.csv("mesocc.pca.leaf.csv", sep=",")
df_group$subs  = factor(df_group$subs, levels=c("Leaf"))
df_group$temp = factor(df_group$temp, levels=c("20","23","26"))
df_group$Day = factor(df_group$Day, levels=c("Day 7","Day 14","Day 21"))
# SEASON SCORE
df_group = subset(df_group, select = c(score, temp, score2, Day) ) # Here you need to input what type of group you want to use

Group_scores <- cbind(Site_tbl, df_group)

# SCORE
Group_scores$score[Group_scores$score=="1"]<-"20" # This is an example, if day is represented by 1 in the df
Group_scores$score[Group_scores$score=="2"]<-"23"
Group_scores$score[Group_scores$score=="3"]<-"26"
Group_scores$score2[Group_scores$score2=="1"]<-"Day 7" # This is an example, if day is represented by 1 in the df
Group_scores$score2[Group_scores$score2=="2"]<-"Day 14"
Group_scores$score2[Group_scores$score2=="3"]<-"Day 21"

# Ordination plot
pca.plot <- ggplot()+
  geom_point(data=Group_scores, aes(x=PC1, y=PC2, fill = temp, shape=Day), size = 4)+
  xlab("PC1 (24% explained var.)")+
  ylab("PC2 (24% explained var.)")+
  scale_fill_manual(values = c("yellow", "orange", "red"))+
  scale_shape_manual(values=c(21,22,23))+
  # geom_segment(data = response_fit, aes(x = 0, xend = PC1, y = 0, yend = PC2), linewidth=1, lineend = "butt",
  #             arrow.fill = "blue", arrow=arrow(length = unit(0.35,"cm"), type = "closed", angle = 20), color="blue") + #for nitrif.md
  #geom_segment(data = response2_fit, aes(x = 0, xend = PC1, y = 0, yend = PC2), linewidth=1, lineend = "butt",
  #             arrow.fill = "green", arrow=arrow(length = unit(0.35,"cm"), type = "closed", angle = 20), color="green") + #for nitrif.dd
  #geom_segment(data = response3_fit, aes(x = 0, xend = PC1, y = 0, yend = PC2), linewidth=1, lineend = "butt",
  #             arrow.fill = "red", arrow=arrow(length = unit(0.35,"cm"), type = "closed", angle = 20), color="red") + #for nitrif.ad
  geom_segment(data=filter(Env_scores, ccatype=="species"), aes(x=0, y=0, xend=PC1, yend=PC2), linewidth=0.4, lineend = "butt", 
               linetype="dashed", arrow.fill = "black", arrow=arrow(length = unit(0.35,"cm"), type = "closed", angle = 20), color="black") +
  #geom_text_repel(aes(x=PC1, y=PC2, label=vgntxt), data=Resp_tbl, seed=400, col="black", size = 4) + # Turn this off if you don't want automatic labels
  coord_fixed(clip = 'off', ylim=c(-2, 2.7), xlim=c(-2, 2)) + # For limiting x and y axes
  #theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #                 panel.border = element_rect(colour = "black", fill=NA, size=0.4))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", size=0.2, fill=NA))+
  theme(axis.ticks=element_line(color="black", linewidth=0.4),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=12, colour="black"), 
        axis.text.x=element_text(size=12), #margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.y=element_text(size=12), #margin = margin(t = 0, r = 10, b = 0, l = 10)), 
        axis.title=element_text(size=12),
        plot.title = element_text(size=12),
        #panel.grid.major = element_blank(), panel.grid.minor  = element_blank(),
        #panel.spacing = unit(2, "lines"),
        #strip.background = element_blank(),
        #strip.text.x = element_text(size=14, colour="black"),
        panel.background = element_rect(fill = NA, color = "black", size=0.4))+#,
  #plot.margin = unit(c(2,1,1,0), "cm"))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(override.aes = list(shape=21)))

pca.plot

# Adonis (or PERMANOVA)
df <-read.csv("mesocc.pca.leaf.csv", sep=",")
df = subset(df , select = c(Day,temp,bank,
                            ER.mgm2h,chla.ugcm2,nh4.bkd,srp.bkd,no3.bkd,doc.mgl,pen.g))

# Betadisperser is a test for multivariate homogeneity of group dispersions (variances), analog to Levene's test
# Good to go if betadispersion between groups are non-significant
dst <- dist(df)
df_bd <- betadisper(dst, df$Day) #p=.11
anova(df_bd)


# Adonis - test froup difference
adonis <- adonis2(dst ~ temp, data = df, method = "euclidian", perm=9999)
adonis
#temp=.4472
#day=.0012

#######Cobble PCA#######
rm(list=ls())
setwd("")
getwd()

list.of.packages<-c("ggplot2", "dplyr", "gtools", "stats", "lme4", "car", "arm",
                    "multcomp","lattice","performance", "MASS", "buildmer", "mgcv", "caret", 
                    "tidyverse", "vegan", "VIM", "reshape2", "ggrepel", "Hmisc")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)>0) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

# Load data
df <- read.csv("mesocc.pca.cob.csv")

df$temp <- factor(df$temp, levels=c("20", "23", "26"))
df$subs  = factor(df$subs, levels=c("Cobble"))

#df$season <- as.character(df$season)

# Subset environmental variables (for PCA)
df_env  = subset(df , select = c(ER.mgm2h,GPP.mgm2h,chla.ugcm2,nh4.bkd,srp.bkd,no3.bkd,doc.mgl))

# Run scaled PCA
df_PCA <- rda(df_env, scale = TRUE)
summary(eigenvals(df_PCA))

# Bind PCs to data
df.corr <- cbind.data.frame(df,
                            scores(df_PCA, display = "sites", scaling = 1))

# pull out correlation variables
data.correlation = subset(df.corr, select = c(Day,Date,Substrate,sub.day,sub.temp,subs,temp,bank,
                                              ER.mgm2h,GPP.mgm2h,chla.ugcm2,nh4.bkd,srp.bkd,no3.bkd,doc.mgl,
                                              nh4.k,srp.k,
                                              PC1,PC2))

correlation <- rcorr(as.matrix(data.correlation[,9:19]), type="pearson")

correlation

# Extract env and response scores and make tables (this don't need any changes)
Resp_sco <- data.frame(summary(df_PCA)$species[,1:3]) # get the species PC1, PC2 & PC3 scores
Site_sco <- data.frame(summary(df_PCA)$sites[,1:3])
Env_sco <- data.frame(summary(df_PCA)$biplot[,1:3])
Resp_tbl <- as_tibble(Resp_sco)
Site_tbl <- as_tibble(Site_sco)
Env_tbl <- as_tibble(Env_sco)
Resp_tbl <- mutate(Resp_tbl, vgntxt=rownames(Resp_sco),
                   ccatype = "species")
Site_tbl <- mutate(Site_tbl, vgntxt=rownames(Site_tbl),
                   ccatype = "sites")
Env_tbl <- mutate(Env_tbl, vgntxt=rownames(Env_sco),
                  ccatype = "bp")

# Use same env variables as in PCA
Env_scores <- Resp_tbl%>% filter(vgntxt %in% c("ER.mgm2h","GPP.mgm2h","chla.ugcm2","nh4.bkd","srp.bkd","no3.bkd","doc.mgl"))

# Use same response variable
Response_scores <- Resp_tbl%>% filter(vgntxt %in% c("srp.k", "nh4.k")) # Same as response variable

# Extract grouping factor i.e., day/night. This is used for coloring points
df_group<-read.csv("mesocc.pca.cob.csv", sep=",")
df_group$subs  = factor(df_group$subs, levels=c("Cobble"))
df_group$temp = factor(df_group$temp, levels=c("20","23","26"))
df_group$Day = factor(df_group$Day, levels=c("Day 7","Day 14","Day 21"))
# SEASON SCORE
df_group = subset(df_group, select = c(score, temp, score2, Day) ) # Here you need to input what type of group you want to use

Group_scores <- cbind(Site_tbl, df_group)

#  SCORE
Group_scores$score[Group_scores$score=="1"]<-"20" # This is an example, if day is represented by 1 in the df
Group_scores$score[Group_scores$score=="2"]<-"23"
Group_scores$score[Group_scores$score=="3"]<-"26"
Group_scores$score2[Group_scores$score2=="1"]<-"Day 7" # This is an example, if day is represented by 1 in the df
Group_scores$score2[Group_scores$score2=="2"]<-"Day 14"
Group_scores$score2[Group_scores$score2=="3"]<-"Day 21"

# Ordination plot
pca.plot <- ggplot()+
  geom_point(data=Group_scores, aes(x=PC1, y=PC2, fill = temp, shape=Day), size = 4)+
  xlab("PC1 (35% explained var.)")+
  ylab("PC2 (25% explained var.)")+
  scale_fill_manual(values = c("yellow", "orange", "red"))+
  scale_shape_manual(values=c(21,22,23))+
  # geom_segment(data = response_fit, aes(x = 0, xend = PC1, y = 0, yend = PC2), linewidth=1, lineend = "butt",
  #             arrow.fill = "blue", arrow=arrow(length = unit(0.35,"cm"), type = "closed", angle = 20), color="blue") + #for nitrif.md
  #geom_segment(data = response2_fit, aes(x = 0, xend = PC1, y = 0, yend = PC2), linewidth=1, lineend = "butt",
  #             arrow.fill = "green", arrow=arrow(length = unit(0.35,"cm"), type = "closed", angle = 20), color="green") + #for nitrif.dd
  #geom_segment(data = response3_fit, aes(x = 0, xend = PC1, y = 0, yend = PC2), linewidth=1, lineend = "butt",
  #             arrow.fill = "red", arrow=arrow(length = unit(0.35,"cm"), type = "closed", angle = 20), color="red") + #for nitrif.ad
  geom_segment(data=filter(Env_scores, ccatype=="species"), aes(x=0, y=0, xend=PC1, yend=PC2), linewidth=0.4, lineend = "butt", 
               linetype="dashed", arrow.fill = "black", arrow=arrow(length = unit(0.35,"cm"), type = "closed", angle = 20), color="black") +
  #geom_text_repel(aes(x=PC1, y=PC2, label=vgntxt), data=Resp_tbl, seed=400, col="black", size = 4) + # Turn this off if you don't want automatic labels
  #stat_ellipse(data=Group_scores, aes(x=PC1, y=PC2, group=Day, color=Day))+
  #scale_color_manual(values = c("#cdf2fa", "#1ec5ea", "#062b34"))+
  coord_fixed(clip = 'off', ylim=c(-2, 2.7), xlim=c(-2, 2)) + # For limiting x and y axes
  #theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #                 panel.border = element_rect(colour = "black", fill=NA, size=0.4))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", size=0.2, fill=NA))+
  theme(axis.ticks=element_line(color="black", linewidth=0.4),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=12, colour="black"), 
        axis.text.x=element_text(size=12), #margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.y=element_text(size=12), #margin = margin(t = 0, r = 10, b = 0, l = 10)), 
        axis.title=element_text(size=12),
        plot.title = element_text(size=12),
        #panel.grid.major = element_blank(), panel.grid.minor  = element_blank(),
        #panel.spacing = unit(2, "lines"),
        #strip.background = element_blank(),
        #strip.text.x = element_text(size=14, colour="black"),
        panel.background = element_rect(fill = NA, color = "black", size=0.4))+#,
  #plot.margin = unit(c(2,1,1,0), "cm"))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(override.aes = list(shape=21)))

pca.plot


# Adonis (or PERMANOVA)
df <-read.csv("mesocc.pca.cob.csv", sep=",")
df = subset(df , select = c(Day,Date,Substrate,sub.day,sub.temp,subs,temp,bank,
                            ER.mgm2h,GPP.mgm2h,chla.ugcm2,nh4.bkd,srp.bkd,no3.bkd,doc.mgl))

# Betadisperser is a test for multivariate homogeneity of group dispersions (variances), analog to Levene's test
# Good to go if betadispersion between groups are non-significant
dst <- dist(df)
df_bd <- betadisper(dst, df$temp)
anova(df_bd)

# Adonis - test group difference
adonis <- adonis2(dst ~ temp, data = df, method = "euclidian", perm=9999)
adonis
#day=.0001
#temp=.6095
######BRMS code######
rm(list=ls())
setwd("")
getwd()
#Packages
library(chron) 
library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggthemes)
library(ggfortify)
library(gridExtra)
library(RColorBrewer)
library(pracma)
library(patchwork)
library(FSA)
library(car)
library(fmsb)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(plyr)
library(rlang)
library(ggpubr)
library(parameters)
library(brms)
library(emmeans)

#dev.off()
#set strings as factors to false
options(stringsAsFactors = FALSE)

#BRMS code from Erik M. Curtis 12/2024
#Updated to work for ANP Meso CC Project

#Use the regression file, remove any regressions with really obvious errors (concentration increased over time, mostly outliers, etc.)
#Run cobble NH4+-N and SRP separate, and leaf NH4+-N and SRP separate. So run this model 4X
#This is an example for cobble NH4

mesodat$temp  = factor(mesodat$temp, levels=c("20", "23","26"))
mesodat$day  = factor(mesodat$day, levels=c("7", "14","21"))

mod <- brm(formula = (ln.nh4) ~ time.1:temp*day + 
             (1 + time.1 | mesocosm:day),     
           
           data = mesodat,
           family = gaussian(),
           prior = c(set_prior("normal(0,10)", class = "b")
           ),
           warmup = 1000, iter = 5000, chains = 6, cores=8,
           control = list(adapt_delta = 0.8), 
           
)

# use emtrends to get estimates + CIs for k (95 here)
Final <- emtrends(mod,c("temp","day"), var = "time.1")
Final

coef(mod)$mesocosm
# export values for k to a data frame so you can plot them
Final.Slopes.CobNH4 <- as.data.frame(Final) # %>% 

Final.Slopes.CobNH4.75 <- Final.Slopes.CobNH4 %>%
  mutate(
    HPD.12.5 = ci(Final, method = "HDI", ci = 0.75)$CI_low,
    HPD.87.5 = ci(Final, method = "HDI", ci = 0.75)$CI_high
  )


######Example GLM Code####
#Example for total chla
glmmTMB_gam <- glmmTMB(Total ~ temp * day + 
                         (1|bank), #Meso.Location (model doesn't converge)
                       family = Gamma(link = "log"), 
                       data = cob.chla, dispformula = ~temp) 

emmmglmm <- emmeans(glmmTMB_gam, ~ temp*day)
plot(emmmglmm)
emmglmm <- emmeans(glmmTMB_gam,list(pairwise~ temp*day), adjust="bonferroni")
emmglmm
simulationOutput <- simulateResiduals(fittedModel = glmmTMB_gam, plot = F)
plot(simulationOutput)