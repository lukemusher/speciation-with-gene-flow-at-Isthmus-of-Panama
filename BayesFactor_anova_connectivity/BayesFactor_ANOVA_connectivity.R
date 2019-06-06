#Author: Luke Musher
setwd("PATH/TO/probabilityNeighborhoods/")

####plot the suitabilities through time####

cur250<-read.csv("250kmProbabilitiesCurrent.csv", header=T)$x
mh250<-read.csv("250kmProbabilitiesMH.csv",header=T)$x
lgm250<-read.csv("250kmProbabilitiesLGM.csv", header=T)$x
lig250<-read.csv("250kmProbabilitiesLIG.csv",header=T)$x

cur50<-read.csv("50kmProbabilitiesCurrent.csv", header=T)$x
mh50<-read.csv("50kmProbabilitiesMH.csv",header=T)$x
lgm50<-read.csv("50kmProbabilitiesLGM.csv", header=T)$x
lig50<-read.csv("50kmProbabilitiesLIG.csv",header=T)$x

par(mfrow=c(1,1))
boxplot(cur50,mh50,lgm50,lig50, horizontal = T,
        col = c("steelblue3", "salmon1","sienna3", "lightblue3"),
        xlab="joint probability",
        names=c("CUR","MH","LGM","LIG" ),
        notch = T, main = "50km neighborhoods")

boxplot(cur250,mh250,lgm250,lig250, horizontal = T,
        col = c("steelblue3", "salmon1","sienna3", "lightblue3"),
        xlab="joint probability",
        names=c("CUR","MH","LGM","LIG" ),
        notch = T, main = "250km neighborhoods")

#################################
###anova with bayes factor
#################################

#install.packages("BayesFactor")
library(BayesFactor)

###START WITH 50KM NEIGHBORHOODS###
data1<-read.csv("../50kmProbNeighb_summarized.csv",header=T)

#First test. Are all categories equal? null hypothesis essentially the same as anova
bf1 = anovaBF(mean.prob~size.time, data=data1)

#next hypothesis, bf2, are lig and cur equal? 
#make new column called 'lig.eg.cur' (LIG equals Current)
data1$lig.eq.cur <- as.character(data1$size.time)

#set LIG and Current probabilities to same time period called LIG/CUR
data1$lig.eq.cur[data1$size.time == "50kmLIG"] <- 'lig/cur'
data1$lig.eq.cur[data1$size.time == "50kmCurrent"] <- 'lig/cur'
data1$lig.eq.cur <- factor(data1$lig.eq.cur)

###Define the new bayesian anova model
bf2 = anovaBF(mean.prob~lig.eq.cur, data = data1)

#next hypothesis, bf3, are lig and MH equal?
data1$lig.eq.MH = as.character(data1$size.time)
data1$lig.eq.MH[data1$size.time == "50kmLIG"] = 'lig/MH'
data1$lig.eq.MH[data1$size.time == "50kmMH"] = 'lig/MH'
data1$lig.eq.MH = factor(data1$lig.eq.MH)

bf3 = anovaBF(mean.prob~lig.eq.MH, data = data1)

#next hypothesis, bf4, are lig and LGM equal?
data1$lig.eq.LGM = as.character(data1$size.time)
data1$lig.eq.LGM[data1$size.time == "50kmLIG"] = 'lig/LGM'
data1$lig.eq.LGM[data1$size.time == "50kmLGM"] = 'lig/LGM'
data1$lig.eq.LGM = factor(data1$lig.eq.LGM)

bf4 = anovaBF(mean.prob~lig.eq.LGM, data = data1)

#next hypothesis, bf5, are MH and LGM equal?
data1$MH.eq.LGM = as.character(data1$size.time)
data1$MH.eq.LGM[data1$size.time == "50kmMH"] = 'MH/LGM'
data1$MH.eq.LGM[data1$size.time == "50kmLGM"] = 'MH/LGM'
data1$MH.eq.LGM = factor(data1$MH.eq.LGM)

bf5 = anovaBF(mean.prob~MH.eq.LGM, data = data1)

#next hypothesis, bf6, are MH and current equal?
data1$MH.eq.cur = as.character(data1$size.time)
data1$MH.eq.cur[data1$size.time == "50kmMH"] = 'MH/CUR'
data1$MH.eq.cur[data1$size.time == "50kmCurrent"] = 'MH/CUR'
data1$MH.eq.cur = factor(data1$MH.eq.cur)

bf6 = anovaBF(mean.prob~MH.eq.cur, data = data1)

###NOW COMPARE ALL MODELS###

bf_all_tests = c(bf1, bf2, bf3, bf4, bf5, bf6)
bf_all_tests

#Bayes factor analysis
#--------------
#[1] size.time  : 5.293502e+160 ±0%
#[2] lig.eq.cur : 3.034144e+63  ±0.03%
#[3] lig.eq.MH  : 1.251744e+159 ±0.01%
#[4] lig.eq.LGM : 1.460272e+156 ±0.01%
#[5] MH.eq.LGM  : 1.399882e+145 ±0.01%
#[6] MH.eq.cur  : 2.104611e+92  ±0.02%

#Against denominator:
#  Intercept only
#---
#  Bayes factor type: BFlinearModel, JZS

###LOOK AT RELATIVE BAYES FACTOR FOR EACH MODEL COMPARED TO MODEL WHERE ALL TIMES VARY###

bf_all_tests[1]/bf_all_tests[4] #all vary vs. LIG == LGM
#[1] size.time : 36250.12 ±0.01%
bf_all_tests[1]/bf_all_tests[5] #all vary vs. LGM == MH
#[1] size.time : 3.781391e+15 ±0.01%
bf_all_tests[1]/bf_all_tests[6] #all vary vs. MH == CUR
#[1] size.time : 2.515193e+68 ±0.02%

###LOOK AT POSTERIOR PROBABILITIES###

#SAMPLE THE POSTERIOR OF THE MODEL
samples = posterior(bf1, iterations = 100000)

## Check order constraint for CUR>LIG>MH>LGM
consistent = (samples[, "size.time-50kmLIG"] > samples[, "size.time-50kmLGM"]) &
  (samples[, "size.time-50kmMH"] > samples[, "size.time-50kmLIG"])&
  (samples[, "size.time-50kmCurrent"] > samples[, "size.time-50kmMH"])
N_consistent = sum(consistent)

bf_restriction_against_full = (N_consistent / 100000) / (1 / 24)
bf_restriction_against_full
#[1] 23.99664 ###VERY SIGNIFICANT
PP = N_consistent/100000
PP
#[1] 0.99986 ###THIS IS A HIGHLY SIGNIFICANT RESULT BASED ON BF AND PP

###Check for shifts in connectiity between time slices:
consistent = (samples[, "size.time-50kmLIG"] > samples[, "size.time-50kmLGM"])
N_consistent = sum(consistent)
PP = (N_consistent / 100000)
PP
#[1] 0.99996

consistent = (samples[, "size.time-50kmMH"] > samples[, "size.time-50kmLGM"])
N_consistent = sum(consistent)
PP = (N_consistent / 100000)
PP
#[1] 0.99997

consistent = (samples[, "size.time-50kmMH"] < samples[, "size.time-50kmCurrent"])
N_consistent = sum(consistent)
PP = (N_consistent / 100000)
PP
#[1] 0.99999

consistent = (samples[, "size.time-50kmLGM"] < samples[, "size.time-50kmCurrent"])
N_consistent = sum(consistent)
PP = (N_consistent / 100000)
PP
#[1] 0.99997

###AS YOU CAN SEE ALL TIME PERIODS VARY WITH SIGNIFICANT SUPPORT###


###NOW LOOK AT 250KM NEIGHBORHOODS###

data2<-read.csv("../250kmProbNeighb_summarized.csv",header=T)

#First test. Are all categories equal? null hypothesis essentially the same as anova
bf1 = anovaBF(mean.prob~size.time, data=data2)
#[1] size.time : 3.674465e+501 ±0% is very small bayes factor and therefore the null hypothesis is rejected

#next hypothesis, bf2, are lig and cur equal?
data2$lig.eq.cur = as.character(data2$size.time)
data2$lig.eq.cur[data2$size.time == "250kmLIG"] = 'lig/cur'
data2$lig.eq.cur[data2$size.time == "250kmCurrent"] = 'lig/cur'
data2$lig.eq.cur = factor(data2$lig.eq.cur)

bf2 = anovaBF(mean.prob~lig.eq.cur, data = data2)

#next hypothesis, bf3, are lig and MH equal?
data2$lig.eq.MH = as.character(data2$size.time)
data2$lig.eq.MH[data2$size.time == "250kmLIG"] = 'lig/MH'
data2$lig.eq.MH[data2$size.time == "250kmMH"] = 'lig/MH'
data2$lig.eq.MH = factor(data2$lig.eq.MH)

bf3 = anovaBF(mean.prob~lig.eq.MH, data = data2)

#next hypothesis, bf4, are lig and LGM equal?
data2$lig.eq.LGM = as.character(data2$size.time)
data2$lig.eq.LGM[data2$size.time == "250kmLIG"] = 'lig/LGM'
data2$lig.eq.LGM[data2$size.time == "250kmLGM"] = 'lig/LGM'
data2$lig.eq.LGM = factor(data2$lig.eq.LGM)

bf4 = anovaBF(mean.prob~lig.eq.LGM, data = data2)

#next hypothesis, bf5, are MH and LGM equal?
data2$MH.eq.LGM = as.character(data2$size.time)
data2$MH.eq.LGM[data2$size.time == "250kmMH"] = 'MH/LGM'
data2$MH.eq.LGM[data2$size.time == "250kmLGM"] = 'MH/LGM'
data2$MH.eq.LGM = factor(data2$MH.eq.LGM)

bf5 = anovaBF(mean.prob~MH.eq.LGM, data = data2)

#next hypothesis, bf6, are MH and current equal?
data2$MH.eq.cur = as.character(data2$size.time)
data2$MH.eq.cur[data2$size.time == "250kmMH"] = 'MH/CUR'
data2$MH.eq.cur[data2$size.time == "250kmCurrent"] = 'MH/CUR'
data2$MH.eq.cur = factor(data2$MH.eq.cur)

bf6 = anovaBF(mean.prob~MH.eq.cur, data = data2)

bf_all_tests = c(bf1, bf2, bf3, bf4, bf5, bf6)
bf_all_tests

# Bayes factor analysis
# --------------
# [1] size.time  : 3.674465e+501 ±0%
# [2] lig.eq.cur : 1.411804e+168 ±0.01%
# [3] lig.eq.MH  : 3.606599e+493 ±0%
# [4] lig.eq.LGM : 1.634135e+486 ±0%
# [5] MH.eq.LGM  : 2.31834e+452  ±0%
# [6] MH.eq.cur  : 2.387684e+250 ±0.01%
#
# Against denominator:
# Intercept only
# ---
# Bayes factor type: BFlinearModel, JZS

bf_all_tests[1]/bf_all_tests[4] #all vary vs. LIG == LGM
#[1] size.time : 2.248568e+15 ±0%
bf_all_tests[1]/bf_all_tests[5] #all vary vs. LGM == MH
#[1] size.time : 1.584955e+49 ±0%
bf_all_tests[1]/bf_all_tests[6] #all vary vs. MH == CUR
#[1] size.time : 1.538924e+251 ±0.01%

samples = posterior(bf1, iterations = 100000)
#head(samples)

## Check order constraint for CUR>LIG>MH>LGM
consistent = (samples[, "size.time-250kmLIG"] > samples[, "size.time-250kmLGM"]) &
  (samples[, "size.time-250kmMH"] > samples[, "size.time-250kmLIG"])&
  (samples[, "size.time-250kmCurrent"] > samples[, "size.time-250kmMH"])
N_consistent = sum(consistent)

bf_restriction_against_full = (N_consistent / 100000) / (1 / 24)
bf_restriction_against_full
#[1] 23.9988
PP = N_consistent/100000
PP
#[1] 0.99995

###Check for shifts in connectiity between time slices:
consistent = (samples[, "size.time-250kmLIG"] > samples[, "size.time-250kmLGM"])
N_consistent = sum(consistent)
PP = (N_consistent / 100000)
PP
#[1] 0.99997

consistent = (samples[, "size.time-250kmMH"] > samples[, "size.time-250kmLGM"])
N_consistent = sum(consistent)
PP = (N_consistent / 100000)
PP
#[1] 0.99996

consistent = (samples[, "size.time-250kmMH"] < samples[, "size.time-250kmCurrent"])
N_consistent = sum(consistent)
PP = (N_consistent / 100000)
PP
#[1] 0.99999

consistent = (samples[, "size.time-250kmLGM"] < samples[, "size.time-250kmCurrent"])
N_consistent = sum(consistent)
PP = (N_consistent / 100000)
PP
#[1] 0.99996
