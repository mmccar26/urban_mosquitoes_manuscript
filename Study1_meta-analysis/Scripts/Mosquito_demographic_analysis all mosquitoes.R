#Yitbarek et al. mosquito demographic analysis
#28 May 2021

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
meta <- file.path(".", "Data", "Mosquito_meta-analysis_data.final.csv")
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
mod0<-rma(yi=LRR,vi=LRR_var, data=mosq.es)
summary(mod0)

#==========Mixed-effects model==================
mod1<-rma.mv(yi=LRR,V=LRR_var, mods=~factor(Taxon)-1, 
             struct="CS",method="REML",digits=4,level = 95, 
             data=mosq.es)
summary(mod1)

#===========Model fit diagnostics===============
par(mfrow=c(1,1)) ## set the plot matrix

##qq norm plots
qqPlot(residuals.rma(mod1), xlab="Theoretical Quantiles", ylab = "Sample Quantiles", line = "quartiles", col.lines = "black", grid = FALSE)

##residual vs fitted plot
plot(fitted.rma(mod1), residuals.rma(mod1), xlab = "Fitted Residuals", ylab = "Residuals") #residuals vs fitted
abline(h=0)

#==========Publication bias analyses====================
#Trim and fill method
res <- rma(LRR, LRR_var, data = mosq.es, method = "FE")
trimfill(res)
funnel(res, legend = TRUE, size =4 )

#Rosenthal's fail-safe number
fsn(LRR, LRR_var, data = mosq.es, type = "Rosenthal")

#=============Forest plots==============================
#By taxon
summary(mod1)
y<-summary(mod1)$b
ci_l<-summary(mod1)$ci.lb
ci_h<-summary(mod1)$ci.ub

msq1<-as.data.frame(cbind(y, ci_l, ci_h))%>%
  add_row(V1 = 0.9155, ci_l = 0.4681, ci_h = 1.3629)%>%
  add_column(Taxon = c("Aedes aegypti", "Aedes albopictus", 
                       "Culex spp.", "Mosquitoes, general", "All mosquitoes"))
  


#ggplot function
msq1%>%
  mutate(Taxon = fct_relevel(Taxon, 
                           "Aedes aegypti", "Aedes albopictus", 
                            "Culex spp.", "Mosquitoes, general", "All mosquitoes"))%>%
  ggplot(aes(x=Taxon, y=V1, ymin=ci_l, ymax=ci_h)) + 
  geom_errorbar(aes(ymin=ci_l, ymax=ci_h), size=1.5, width=0.1)+
  geom_point(size = 4)+
  geom_pointrange()+
  coord_flip(ylim=c(-2, 2))+
  scale_y_continuous(breaks = seq(-2, 2, by = 1))+
  xlab(NULL) +
  ylab(NULL)+
  geom_hline(aes(yintercept=0))+
  theme(axis.text = element_text(size=12, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))

##all added
mosq.es%>%
ggplot(aes(x=Taxon, y=LRR)) + 
  geom_point(alpha = .4)+
  geom_point(data = msq1$v1, size = 2)+
  #geom_pointrange()+
  #coord_flip(ylim=c(-3, 3))+
  xlab(NULL) +
  ylab(NULL)+
  geom_hline(aes(yintercept=0))+
  theme(axis.text = element_text(size=12, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))



#all mosquitoes
#By taxon
summary(mod0)
y<-summary(mod0)$b
ci_l<-summary(mod0)$ci.lb
ci_h<-summary(mod0)$ci.ub

msq2<-as.data.frame(cbind(y, ci_l, ci_h))%>%
  add_column(Taxon = c("All Mosquitoes"))

#ggplot function
msq2%>%
  ggplot(aes(x= Taxon, y = V1, ymin=ci_l, ymax=ci_h)) + 
  geom_errorbar(aes(ymin=ci_l, ymax=ci_h), size=1.5, width=0.1)+
  geom_point(shape=19, size=4)+
  geom_pointrange()+
  coord_flip(ylim=c(-2, 2))+
  xlab(NULL) +
  ylab(NULL)+
  geom_hline(aes(yintercept=0))+
  theme(axis.text = element_text(size=12, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))


