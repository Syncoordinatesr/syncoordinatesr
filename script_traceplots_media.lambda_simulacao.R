#Traceplots e media.lambda Simulação

#set.seed()
library(rgdal)
library(ggplot2)
library(rgeos)
library(plyr)
library(dplyr)
library(scales)
library(epiDisplay)
library(patchwork)
library(sf)
library(FNN)
library(ars)
library(fields)
library(mvtnorm)
library(spdep)
library(spBayes)
library(class)
library(HDInterval)
library(coda)
library(MfUSampler)
#install.packages("devtools")
library(devtools)
#devtools::install_github("leogalhardo/syncoordinatesr")
library(syncoordinatesr)
library(yarrr)
library(ggpubr)
library(gridExtra)

#####
#FUNÇÃO: prepare_data
prepara_dados = prepare_data(dados, coords, limits = c(0,10,0,10),
                             grid = 10, continuous = 4)

#####
#FUNÇÃO: syn_mcmc
mcmc_dados = syn_mcmc(dados, coords, limits = c(0,10,0,10),
                      grid = 10, S = 10000, burn = 1000,
                      continuous = 4, 
                      spatial_beta = FALSE,
                      return_parameters = TRUE)

#####
#Parâmetros do mcmc
grid = 10
S = mcmc_dados[["S"]]
burn = mcmc_dados[["burn"]]
lambda = mcmc_dados[["lambda"]]
media.lambda = mcmc_dados[["media.lambda"]]
alfa = mcmc_dados[["alfa"]]
mu = mcmc_dados[["mu"]]
theta = mcmc_dados[["theta"]]
tau.theta = mcmc_dados[["tau.theta"]]
phi = mcmc_dados[["phi"]]
tau.phi = mcmc_dados[["tau.phi"]]
epsilon = mcmc_dados[["epsilon"]]
tau.e = mcmc_dados[["tau.e"]]


#####
#Gráficos de convergência

#traceplots dos parâmetros do modelo
pdf("convergência.pdf", width=20, height=35, pointsize=30)
par(mfrow=c(4,3))
plot(mu[1:S], type = "l", xlab = "Iterações", ylab = "mu")
plot(alfa[1,1:S], type = "l", xlab = "Iterações", ylab = "alfa (1)")
plot(alfa[2,1:S], type = "l", xlab = "Iterações", ylab = "alfa (2)")
plot(alfa[3,1:S], type = "l", xlab = "Iterações", ylab = "alfa (3)")
plot(alfa[4,1:S], type = "l", xlab = "Iterações", ylab = "alfa (4)")
plot(tau.theta[1:S], type = "l", xlab = "Iterações", ylab = "tau.theta")
plot(tau.phi[1,1:S], type = "l", xlab = "Iterações", ylab = "tau.phi (1)")
plot(tau.phi[2,1:S], type = "l", xlab = "Iterações", ylab = "tau.phi (2)")
plot(tau.phi[3,1:S], type = "l", xlab = "Iterações", ylab = "tau.phi (3)")
plot(tau.phi[4,1:S], type = "l", xlab = "Iterações", ylab = "tau.phi (4)")
plot(tau.e[1:S], type = "l", xlab = "Iterações", ylab = "tau.e")
dev.off()

#traceplots dos parâmetros do modelo
pdf("convergência_sem_burnin.pdf", width=20, height=35, pointsize=30)
par(mfrow=c(4,3))
plot(mu[1:S], type = "l", xlab = "Iterações", ylab = "mu")
plot(alfa[1,1001:S], type = "l", xlab = "Iterações", ylab = "alfa (1)")
plot(alfa[2,1001:S], type = "l", xlab = "Iterações", ylab = "alfa (2)")
plot(alfa[3,1001:S], type = "l", xlab = "Iterações", ylab = "alfa (3)")
plot(alfa[4,1001:S], type = "l", xlab = "Iterações", ylab = "alfa (4)")
plot(tau.theta[1001:S], type = "l", xlab = "Iterações", ylab = "tau.theta")
plot(tau.phi[1,1001:S], type = "l", xlab = "Iterações", ylab = "tau.phi (1)")
plot(tau.phi[2,1001:S], type = "l", xlab = "Iterações", ylab = "tau.phi (2)")
plot(tau.phi[3,1001:S], type = "l", xlab = "Iterações", ylab = "tau.phi (3)")
plot(tau.phi[4,1001:S], type = "l", xlab = "Iterações", ylab = "tau.phi (4)")
plot(tau.e[1001:S], type = "l", xlab = "Iterações", ylab = "tau.e")
dev.off()

#plot de lambda(56)
pdf("lambda(56)_sem_burnin.pdf", width=20, height=35, pointsize=30)
par(mfrow=c(4,3))
plot(lambda[56,1001:S,1], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 1")
plot(lambda[56,1001:S,2], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 2")
plot(lambda[56,1001:S,3], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 3")
plot(lambda[56,1001:S,4], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 4")
plot(lambda[56,1001:S,5], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 5")
plot(lambda[56,1001:S,6], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 6")
plot(lambda[56,1001:S,7], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 7")
plot(lambda[56,1001:S,8], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 8")
plot(lambda[56,1001:S,9], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 9")
plot(lambda[56,1001:S,10], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 10")
plot(lambda[56,1001:S,11], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 11")
plot(lambda[56,1001:S,12], type = "l", xlab = "Iterações", ylab = "Lambda (56) combinação 12")
dev.off()


#####
#plot da avaliação do modelo
corDegrade <- colorRampPalette(c("light yellow","red"))

pdf("media.lambda.pdf", width=30, height=40, pointsize=40)
par(mfrow=c(4,3))

image(lonvec,latvec,matrix(media.lambda[,1],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=0,x1=1,x2=1"))
points(dados$lon[comb==1],dados$lat[comb==1])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,2],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=0,x1=1,x2=2"))
points(dados$lon[comb==2],dados$lat[comb==2])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,3],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=0,x1=1,x2=3"))
points(dados$lon[comb==3],dados$lat[comb==3])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,4],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=0,x1=2,x2=1"))
points(dados$lon[comb==4],dados$lat[comb==4])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,5],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=0,x1=2,x2=2"))
points(dados$lon[comb==5],dados$lat[comb==5])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,6],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=0,x1=2,x2=3"))
points(dados$lon[comb==6],dados$lat[comb==6])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,7],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=1,x1=1,x2=1"))
points(dados$lon[comb==7],dados$lat[comb==7])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,8],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=1,x1=1,x2=2"))
points(dados$lon[comb==8],dados$lat[comb==8])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,9],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=1,x1=1,x2=3"))
points(dados$lon[comb==9],dados$lat[comb==9])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,10],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=1,x1=2,x2=1"))
points(dados$lon[comb==10],dados$lat[comb==10])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,11],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=1,x1=2,x2=2"))
points(dados$lon[comb==11],dados$lat[comb==11])
rect(4, 4, 6, 5, border = "black")

image(lonvec,latvec,matrix(media.lambda[,12],nrow=grid,ncol=grid),col=corDegrade(10),main=paste("Intensidade y=1,x1=2,x2=3"))
points(dados$lon[comb==12],dados$lat[comb==12])
rect(4, 4, 6, 5, border = "black")

dev.off()
