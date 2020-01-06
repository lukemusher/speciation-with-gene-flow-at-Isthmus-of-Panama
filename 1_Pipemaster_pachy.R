###############################################################################
####Code written by Greg Thom, modified from code written by Marcelo Gehara####
###############################################################################

library(PipeMaster)
library(ggplot2)
library(caret)
library(doMC)
library(ggpubr)
setwd("~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/")

###Drawing the models
m1 <- main.menu(m1)
m2 <- main.menu(m2)
m3 <- main.menu(m3)
m4 <- main.menu(m4)


pop.assign <- read.delim("pachy_pops.txt", header = FALSE)


m1 <- get.data.structure(m1,path.to.fasta = "~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/fasta_trim_Pipe_master/", pop.assign=pop.assign)

##observed_sum_stats
stat_pachy <- obs.sumstat.ngs(model = m1, path.to.fasta = "~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/fasta_trim_Pipe_master/", pop.assign = pop.assign, moments = F)

###Running the models###

m2 <- get.data.structure(m2,path.to.fasta = "~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/fasta_trim_Pipe_master/", pop.assign=pop.assign)
m3 <- get.data.structure(m3,path.to.fasta = "~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/fasta_trim_Pipe_master/", pop.assign=pop.assign)
m4 <- get.data.structure(m4,path.to.fasta = "~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/fasta_trim_Pipe_master/", pop.assign=pop.assign)

###m1###  
sim.msABC.sumstat(m1, path="~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/",
                  #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                  nsim.blocks = 10, use.alpha = F,
                  output.name = "m1_pachy",
                  append.sims = F, ncores=8, block.size = 300)
###m2###
sim.msABC.sumstat(m2, path="~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/",
                  #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                  nsim.blocks = 10, use.alpha = F,
                  output.name = "m2_pachy",
                  append.sims = F, ncores=8, block.size = 300)

###m3###
sim.msABC.sumstat(m3, path="~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/",
                  #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                  nsim.blocks = 10, use.alpha = F,
                  output.name = "m3_pachy",
                  append.sims = F, ncores=5, block.size = 480)

###m4###
sim.msABC.sumstat(m4, path="~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/",
                  #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                  nsim.blocks = 10, use.alpha = F,
                  output.name = "m4_pachy",
                  append.sims = F, ncores=8, block.size = 300)

###########################PCA
##clustering simulations
setwd("~/Dropbox (ecoevounifesp)/Publicacoes em andamento/pachyramphus _luke/Pipe_master_pachy/Pipe_master_pachy_2/")
m1.sim <- read.delim("SIMS_m1_pachy.txt")
m2.sim <- read.delim("SIMS_m2_pachy.txt")
m3.sim <- read.delim("SIMS_m3_pachy.txt")
m4.sim <- read.delim("SIMS_m4_pachy.txt")


data <- c(rep("m1",nrow(m1.sim)),
          rep("m2",nrow(m2.sim)),
          rep("m3",nrow(m3.sim)),
          rep("m4",nrow(m4.sim)))
          
          

data_pca <- c(rep("m1",1000),
          rep("m2",1000),
          rep("m3",1000),
          rep("m4",1000))
          


models <- rbind(m1.sim[,6:ncol(m1.sim)],
                m2.sim[,8:ncol(m2.sim)],
                m3.sim[,12:ncol(m3.sim)],
                m4.sim[,12:ncol(m4.sim)])
                

models_pca <- rbind(m1.sim[1:1000,6:ncol(m1.sim)],
                m2.sim[1:1000,8:ncol(m2.sim)],
                m3.sim[1:1000,12:ncol(m3.sim)],
                m4.sim[1:1000,12:ncol(m4.sim)])
               


  models_pca <- models_pca[,-grep("Fay",colnames(models_pca))]
  
  models_pca <- models_pca[,-grep("thomson",colnames(models_pca))]
  
  
 
  data_pca<-data_pca[complete.cases(models_pca)]
  models_pca<-models_pca[complete.cases(models_pca),]
  
  
  models <- models[,-grep("Fay",colnames(models))]
  
  models <- models[,-grep("thomson",colnames(models))]
  
  
stat_pachy<-stat_pachy[colnames(stat_pachy) %in% colnames(models_pca)]

  PCA <- prcomp(rbind(models_pca,stat_pachy), center = T, scale. = T, retx=T)
  scores <- data.frame(PCA$x[,1:ncol(PCA$x)])
  
p1x2<-ggplot(scores, aes(x=PC1, y=PC2))+ theme(legend.position = "none") +  theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(2,2,2,2,4))+
  scale_colour_manual(values = c("yellow", "green", "blue", "red", "black")) + scale_shape_manual(values = c(1,2,3,4,16))

p1x3<-ggplot(scores, aes(x=PC1, y=PC3))+ theme(legend.position = "none") +  theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(2,2,2,2,4))+
  scale_colour_manual(values = c("yellow", "green", "blue", "red", "black")) + scale_shape_manual(values = c(1,2,3,4,16))

p2x3<-ggplot(scores, aes(x=PC2, y=PC3))+ theme(legend.position = "none") +  theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(2,2,2,2,4))+
  scale_colour_manual(values = c("yellow", "green", "blue", "red", "black")) + scale_shape_manual(values = c(1,2,3,4,16))

PCA_Plot <- ggarrange(p1x2,p1x3,p2x3, legend = "left", common.legend = T, ncol = 3, nrow = 1)
pdf("Pachyramphus.pdf", width=12, height=3)
PCA_Plot
dev.off()

##ABC
library(abc)
target= stat_pachy

a_0.1 <- postpr(target=target, index=data, 
            sumstat= models, tol=0.1, method="neuralnet", kernel="epanechnikov",
            sizenet = 20, maxit = 1000,MaxNWts = 15000, hcorr = TRUE , corr=FALSE, trace = TRUE)

a_0.05 <- postpr(target=target, index=data, 
                sumstat= models, tol=0.05, method="neuralnet", kernel="epanechnikov",
                sizenet = 20, maxit = 1000,MaxNWts = 15000, hcorr = TRUE , corr=FALSE, trace = TRUE)

a_0.01 <- postpr(target=target, index=data, 
                sumstat= models, tol=0.01, method="neuralnet", kernel="epanechnikov",
                sizenet = 20, maxit = 1000,MaxNWts = 15000, hcorr = TRUE , corr=FALSE, trace = TRUE)

summary(a_0.1)
summary(a_0.05)
summary(a_0.01)


sumstat <- m4.sim[,colnames(m4.sim) %in% colnames(models)]

x0.1 <- abc(stat_pachy, param = m4.sim[,1:11], sumstat = sumstat, tol=0.1, sizenet = 20, maxit = 1000, MaxNWts = 15000, method ="neuralnet")
x0.05 <- abc(stat_pachy, param = m4.sim[,1:11], sumstat = sumstat, tol=0.05, sizenet = 20, maxit = 1000, MaxNWts = 15000, method ="neuralnet")
x0.01 <- abc(stat_pachy, param = m4.sim[,1:11], sumstat = sumstat, tol=0.01, sizenet = 20, maxit = 1000, MaxNWts = 15000, method ="neuralnet")


t0.1 <-summary(x0.1)
t0.05 <-summary(x0.05)
t0.01 <-summary(x0.01)


###########CV4###############
#arquivo com os parametros
library(plyr)

sim_param <- rbind.fill(m1.sim[,1:5],
                        m2.sim[,1:7],
                        m3.sim[,1:11],
                        m4.sim[,1:11])

sim_param[is.na(sim_param)] <- 0

cv4postpr_neuralnet_0.1=cv4postpr(param=sim_param, sumstat=models, nval=100, tol=0.1, method="neuralnet", kernel="epanechnikov",
                                   sizenet = 20, maxit = 1000, MaxNWts = 15000, numnet = 10, index=data, statistic='mean')
z_0.1 <- summary(cv4postpr_neuralnet_0.1)

cv4postpr_neuralnet_0.05=cv4postpr(param=sim_param, sumstat=models, nval=100, tol=0.05, method="neuralnet", kernel="epanechnikov",
                              sizenet = 20, maxit = 1000, MaxNWts = 15000, numnet = 10, index=data, statistic='mean')
z_0.05 <- summary(cv4postpr_neuralnet_0.05)

cv4postpr_neuralnet_0.01=cv4postpr(param=sim_param, sumstat=models, nval=100, tol=0.01, method="neuralnet", kernel="epanechnikov",
                              sizenet = 20, maxit = 1000, MaxNWts = 15000, numnet = 10, index=data, statistic='mean')
z_0.01 <- summary(cv4postpr_neuralnet_0.01)


library(lattice)
library(RColorBrewer)
breaks <- seq(0, 1, by=0.05)
breaks2 <- seq(-1, 100, by=5.0)
myColours2 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))


data0.1 <- round(t(as.matrix(z_0.1$probs[[1]])), 2)
myPanel0.1 <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, data0.1[cbind(x,y)]) ## use handy matrix indexing
}
data2_0.1 <- t(as.matrix(z_0.1$conf.matrix[[1]]))
conf0.1 <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, data2_0.1[cbind(x,y)]) ## use handy matrix indexing
}

data0.05 <- round(t(as.matrix(z_0.05$probs[[1]])), 2)
myPanel0.05 <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, data0.05[cbind(x,y)]) ## use handy matrix indexing
}
data2_0.05 <- t(as.matrix(z_0.05$conf.matrix[[1]]))
conf0.05 <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, data2_0.05[cbind(x,y)]) ## use handy matrix indexing
}

data0.01 <- round(t(as.matrix(z_0.01$probs[[1]])), 2)
myPanel0.01 <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, data0.01[cbind(x,y)]) ## use handy matrix indexing
}
data2_0.01 <- t(as.matrix(z_0.01$conf.matrix[[1]]))
conf0.01 <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, data2_0.01[cbind(x,y)]) ## use handy matrix indexing
}

CV_confussion0.1 <- levelplot(t(as.matrix(z_0.1$probs[[1]])), col.regions = myColours2, at=breaks, panel=myPanel0.1, colorkey=list(labels=list(cex=1)), scales=list(x=list(cex=1),y=list(cex=1)), xlab=list("Simulated", cex=1.2), ylab=list("Observed", cex=1.2), cex.lab=1.5)
CV_confussion0.05 <- levelplot(t(as.matrix(z_0.05$probs[[1]])), col.regions = myColours2, at=breaks, panel=myPanel0.05, colorkey=list(labels=list(cex=1)), scales=list(x=list(cex=1),y=list(cex=1)), xlab=list("Simulated", cex=1.2), ylab=list("Observed", cex=1.2), cex.lab=1.5)
CV_confussion0.01 <- levelplot(t(as.matrix(z_0.01$probs[[1]])), col.regions = myColours2, at=breaks, panel=myPanel0.01, colorkey=list(labels=list(cex=1)), scales=list(x=list(cex=1),y=list(cex=1)), xlab=list("Simulated", cex=1.2), ylab=list("Observed", cex=1.2), cex.lab=1.5)
CV_probs0.1 <- levelplot(t(as.matrix(z_0.1$conf.matrix[[1]])), col.regions = myColours2, at=breaks2, panel=conf0.1, colorkey=list(labels=list(cex=1)), scales=list(x=list(cex=1),y=list(cex=1)), xlab=list("Simulated", cex=1.2), ylab=list("Observed", cex=1.2), cex.lab=1.5)
CV_probs0.05 <- levelplot(t(as.matrix(z_0.05$conf.matrix[[1]])), col.regions = myColours2, at=breaks2, panel=conf0.05, colorkey=list(labels=list(cex=1)), scales=list(x=list(cex=1),y=list(cex=1)), xlab=list("Simulated", cex=1.2), ylab=list("Observed", cex=1.2), cex.lab=1.5)
CV_probs0.01 <- levelplot(t(as.matrix(z_0.01$conf.matrix[[1]])), col.regions = myColours2, at=breaks2, panel=conf0.01, colorkey=list(labels=list(cex=1)), scales=list(x=list(cex=1),y=list(cex=1)), xlab=list("Simulated", cex=1.2), ylab=list("Observed", cex=1.2), cex.lab=1.5)

library(gridExtra)
pdf("pachy_heatmap.pdf", width=12, height=6)
grid.arrange(CV_confussion0.1,CV_confussion0.05,CV_confussion0.01,CV_probs0.1,CV_probs0.05,CV_probs0.01, ncol=3)
dev.off()

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

####Cross_validation for parameters######

cv4abc.object=cv4abc(param=m4.sim[,1:11], sumstat=sumstat, nval=100, tol=c(0.1, 0.05, 0.01), method='neuralnet', 
                     sizenet = 20, maxit = 2000, MaxNWts = 10000, numnet = 10, index=data, statistic='mean')

summary(cv4abc.object)

correlation=cor(cv4abc.object$true,cv4abc.object$estim$tol0.05)


e0.05 <- cv4abc.object$estim$tol0.05
colnames(e0.05) <- paste("sim", colnames(e0.05), sep = "")
e0.05 <- cbind(cv4abc.object$true, e0.05)
colnames(cv4abc.object$true)

Ne0.pop1 <- ggscatter(e0.05, x = "Ne0.pop1", y = "simNe0.pop1", title = "Ne0.pop1", size = 2, xlab ="Simulated", ylab = "Observed",
          add = "reg.line",                                
          conf.int = TRUE,                              
          add.params = list(color = "blue", fill = "lightgray")
)+
  stat_cor(method = "pearson", cex=5)  + font("xlab", size = 0)+
  font("ylab", size = 0) +  font("title", size = 24) +  font("xy.text", size = 15)

Ne0.pop2 <- ggscatter(e0.05, x = "Ne0.pop2", y = "simNe0.pop2", title = "Ne0.pop2", size = 2, xlab ="Simulated", ylab = "Observed",
                      add = "reg.line",                                
                      conf.int = TRUE,                              
                      add.params = list(color = "blue", fill = "lightgray")
)+
  stat_cor(method = "pearson", cex=5)  + font("xlab", size = 0)+
  font("ylab", size = 0) +  font("title", size = 24) +  font("xy.text", size = 15)

join1_2 <- ggscatter(e0.05, x = "join1_2", y = "simjoin1_2", title = "join1_2", size = 2, xlab ="Simulated", ylab = "Observed",
                      add = "reg.line",                                
                      conf.int = TRUE,                              
                      add.params = list(color = "blue", fill = "lightgray")
)+
  stat_cor(method = "pearson", cex=5)  + font("xlab", size = 0)+
  font("ylab", size = 0) +  font("title", size = 24) +  font("xy.text", size = 15)

mig0.1_2 <- ggscatter(e0.05, x = "mig0.1_2", y = "simmig0.1_2", title = "mig0.1_2", size = 2, xlab ="Simulated", ylab = "Observed",
                      add = "reg.line",                                
                      conf.int = TRUE,                              
                      add.params = list(color = "blue", fill = "lightgray")
)+
  stat_cor(method = "pearson", cex=5)  + font("xlab", size = 0)+
  font("ylab", size = 0) +  font("title", size = 24) +  font("xy.text", size = 15)

mig0.2_1 <- ggscatter(e0.05, x = "mig0.2_1", y = "simmig0.2_1", title = "mig0.2_1", size = 2, xlab ="Simulated", ylab = "Observed",
                      add = "reg.line",                                
                      conf.int = TRUE,                              
                      add.params = list(color = "blue", fill = "lightgray")
)+
  stat_cor(method = "pearson", cex=5)  + font("xlab", size = 0)+
  font("ylab", size = 0) +  font("title", size = 24) +  font("xy.text", size = 15)

mean.rate <- ggscatter(e0.05, x = "mean.rate", y = "simmean.rate", title = "mean.rate", size = 2, xlab ="Simulated", ylab = "Observed",
                      add = "reg.line",                                
                      conf.int = TRUE,                              
                      add.params = list(color = "blue", fill = "lightgray")
)+
  stat_cor(method = "pearson", cex=5)  + font("xlab", size = 0)+
  font("ylab", size = 0) +  font("title", size = 24) +  font("xy.text", size = 15)

sd.rate <- ggscatter(e0.05, x = "sd.rate", y = "simsd.rate", title = "sd.rate", size = 2, xlab ="Simulated", ylab = "Observed",
                      add = "reg.line",                                
                      conf.int = TRUE,                              
                      add.params = list(color = "blue", fill = "lightgray")
)+
  stat_cor(method = "pearson", cex=5)  + font("xlab", size = 0)+
  font("ylab", size = 0) +  font("title", size = 24) +  font("xy.text", size = 15)

library(ggcorrplot)
correlogram <- ggcorrplot(correlation, outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"),
           lab = T)

pdf("Correlation.pdf", width = 18, height = 25)
Correlation <- grid.arrange(Ne0.pop1, Ne0.pop2, join1_2, mig0.1_2, mig0.2_1, mean.rate, sd.rate, ncol = 3, nrow = 3)
dev.off()

pdf("correlogram.pdf", width = 5, height = 5)
correlogram
dev.off()

