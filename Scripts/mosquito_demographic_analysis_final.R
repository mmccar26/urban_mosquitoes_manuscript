#Yitbarek et al. mosquito demographic analysis
#15 December 2022

#=====load libraries======
library(ggplot2)
library(metafor)
library(plyr)
library(lme4)
library(car)
library(dplyr)
library(tidyverse)
library(broom)
library(Rmisc)

#=====import data======
#relative pathname
meta <- file.path(".", "Data", "urban_mosquito_analysis_final.csv")
print(meta)

#import data
mosquito <- read_csv(meta)

#Calculate log Response Ratio
mosq.es<-escalc(measure="ROM",
                  m1i= Mean.low, m2i=Mean.high,
                  sd1i=Std.dev.low, sd2i=Std.dev.high,
                  n1i=SS.low, n2i= SS.high, 
                  data=mosquito,
                  var.names=c("LRR","LRR_var"),digits=4)

##Number of observations per Taxon
mosq.es%>%
  group_by(Taxon)%>%
  tally()

#==========Random-effects model=================
mod0<-rma.mv(yi=LRR,V=LRR_var, random = list(~1|Study, ~1|ID), data=mosq.es)
summary(mod0)
anova(mod0)

#===========Model fit diagnostics===============
par(mfrow=c(1,1)) ## set the plot matrix

##qq norm plots
qqPlot(residuals.rma(mod0), xlab="Theoretical Quantiles", ylab = "Sample Quantiles", line = "quartiles", col.lines = "black", grid = FALSE)

#==========Publication bias analyses====================
#Trim and fill method
res <- rma(LRR, LRR_var, data = mosq.es, method = "FE")
trimfill(res)
funnel(res, legend = TRUE, size = 4)

#Rosenthal's fail-safe number
fsn(LRR, LRR_var, data = mosq.es, type = "Rosenthal")

#==========Mixed-effects model==================

#remove Ae. japonicus due to low sample size
mosq.es<-
  mosq.es%>%
  filter(Taxon != "Aedes japonicus")

#run separate mixed-effects models
mod1<-rma.mv(yi=LRR,V=LRR_var, mods=~factor(Taxon)-1, 
             struct="CS",method="REML",digits=4,level = 95, 
             data=mosq.es)
summary(mod1)


mod2<-rma.mv(yi=LRR,V=LRR_var, mods=~factor(Taxon)-1, random = list(~1|Study),
                 struct="CS",method="REML",digits=4,level = 95, 
                 data=mosq.es)
summary(mod2)

mod3<-rma.mv(yi=LRR,V=LRR_var, mods=~factor(Taxon)-1, random = list(~1|ID),
             struct="CS",method="REML",digits=4,level = 95, 
             data=mosq.es)
summary(mod3)

mod4<-rma.mv(yi=LRR,V=LRR_var, mods=~factor(Taxon)-1, random = list(~1|Study, ~1|ID),
             struct="CS",method="REML",digits=4,level = 95, 
             data=mosq.es)
summary(mod4)


#=============Forest plots==============================
#By taxon
summary(mod2)
y<-summary(mod2)$b
ci_l<-summary(mod2)$ci.lb
ci_h<-summary(mod2)$ci.ub

msq1<-as.data.frame(cbind(y, ci_l, ci_h))%>%
      add_row(V1 =  0.4881, ci_l = 0.0081, ci_h = 0.9681)%>%
add_column(Taxon = c("Aedes aegypti", "Aedes albopictus", 
                       "Culex spp.", "Mosquitoes, non-dist.", "All mosquitoes"))
  
  #ggplot function
msq1%>%
  mutate(Taxon = fct_relevel(Taxon, 
                           "Mosquitoes, non-dist.","Culex spp.","Aedes albopictus", "Aedes aegypti", "All mosquitoes"))%>%
  ggplot(aes(x=Taxon, y=V1, ymin=ci_l, ymax=ci_h, shape = Taxon, color = Taxon)) + 
  geom_errorbar(aes(ymin=ci_l, ymax=ci_h), size=1.5, width=0.1)+
  geom_point(size = 5)+
  scale_shape_manual(values = c(18, 18, 18, 18, 16)) +
  scale_color_manual(values = c("black", "black", "black", "black", "red"))+
  geom_pointrange()+
  coord_flip(ylim=c(-3, 3))+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  xlab(NULL) +
  ylab(NULL)+
  geom_hline(aes(yintercept=0))+
  theme(axis.text = element_text(size=16, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))

