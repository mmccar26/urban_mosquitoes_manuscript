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
meta <- file.path(".", "Data", "Mosquito_meta-analysis_covariate_final.csv")
print(meta)

#import data
cov <- read_csv(meta)

#Calculate log Response Ratio
covar<-escalc(measure="ROM",
                m1i= Mean.low, m2i=Mean.high,
                sd1i=Std.dev.low, sd2i=Std.dev.high,
                n1i=SS.low, n2i= SS.high, 
                data=cov,
                var.names=c("LRR","LRR_var"),digits=4)

##Number of observations per covariate
covar%>%
  group_by(Covariate)%>%
  tally()

#==========Publication bias analyses====================
#Trim and fill method
res <- rma(LRR, LRR_var, data = covar)
trimfill(res)
funnel(res, legend = TRUE, size = 4)

#Rosenthal's fail-safe number
fsn(LRR, LRR_var, data = covar, type = "Rosenthal")

#==========Mixed-effects model==================
mod1<-rma.mv(yi=LRR,V=LRR_var, mods=~factor(Covariate)-1, random = list (~1|Study, ~1|ID),
             struct="CS", method="REML",digits=4,level = 95, 
             data=covar)
summary(mod1)

mod2<-rma.mv(yi=LRR,V=LRR_var, mods=~factor(Covariate)-1, random = list (~1|Study),
             struct="CS", method="REML",digits=4,level = 95, 
             data=covar)
summary(mod2)

mod3<-rma.mv(yi=LRR,V=LRR_var, mods=~factor(Covariate)-1, random = list (~1|ID),
             struct="CS", method="REML",digits=4,level = 95, 
             data=covar)
summary(mod3)

mod4<-rma.mv(yi=LRR,V=LRR_var, mods=~factor(Covariate)-1,
             struct="CS", method="REML",digits=4,level = 95, 
             data=covar)
summary(mod4)


#=============Forest plots==============================
#By covariate
summary(mod3)
y<-summary(mod3)$b
ci_l<-summary(mod3)$ci.lb
ci_h<-summary(mod3)$ci.ub

msq1<-as.data.frame(cbind(y, ci_l, ci_h))%>%
  add_column(Covariate = c("Abandoned buildings", "Containers", 
                       "Education", "Vegetation"))

#ggplot function
msq1%>%
  mutate(Covariate = fct_relevel(Covariate, 
                             "Vegetation","Education",  
                              "Containers","Abandoned buildings"))%>%
  ggplot(aes(x=Covariate, y=V1, ymin=ci_l, ymax=ci_h, shape = Covariate, color = Covariate)) + 
  geom_errorbar(aes(ymin=ci_l, ymax=ci_h), size=1.5, width=0.1)+
  geom_point(size = 5)+
  scale_shape_manual(values = c(16, 16, 16, 16)) +
  scale_color_manual(values = c("black", "black", "black", "black"))+
  geom_pointrange()+
  coord_flip(ylim=c(-2, 2))+
  scale_y_continuous(breaks = seq(-2, 2, by = 1))+
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
