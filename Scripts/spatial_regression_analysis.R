#Yitbarek et al. mosquito spatial analysis
#September 2021

#=====load libraries===================
library(spBayes)
library(mapproj)
library(maps)
library(geoR)

#=====import data======================

# Load the mosquito data
data(mosquito)

#=====Longitude and Latitude============

# longitude and latitude coordinates and years.
coord = cbind(mosquito$long,mosquito$lat)
dates = mosquito$year

#=====Gaussian process gradient model====

#MCMC tuning parameters (default values)
amcmc=list(n.batch=10,batch.length=100,accept.rate=0.3)

#=====Bayesian spatial model=============
attach(localgrad) #localgrad function

#Fit Bayesian spatial regression model
out.grad = localgrad(dates=dates,coord.longlat=coord, Albers=c(29.5,45.5), n.samp=1000, amcmc=amcmc)
class(out.grad)

spatial <- spLM(year~mhi + ttrees + zmean + zsd + zabandon + cii + cip + ciib, n.samples = 10000)

#Export data to csv file
write.csv(spatial$p.beta.tauSq.samples, "/Users/senayyitbarek/Desktop/bayesianresults.csv")

#Credibility intervals: 95 %
library(condo)
bayes = read.csv(file = 'data/bayesianresults.csv')
HPDinterval(as.mcmc(bayes$Xmhi))


#=====Summary and visualization==========
summary(out.grad) 

#plotgrad functions to draw vector field plot
plotgrad(out.grad,cex=1,pch=".",database="state")





