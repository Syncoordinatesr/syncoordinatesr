

#carregamento dos pacotes
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
library(readr)


#carregamento do banco de dados para simulação
dados_simulacao <- read_csv("dados_originais_completo_simulações.csv")
dados_sem_coord <-dados_simulacao[,2:5]
coordenadas <- dados_simulacao[,c(6,7)]

#carregamento utills_functions

# Chame a função prepare_data
prepara_dados <- prepare_data(
  dataset = dados_sem_coord,
  grid = 10,
  continuous = 4,
  limits = c(0, 10, 0, 10),
  coord = coordenadas
)


dados_sem_coord<- as.data.frame(dados_sem_coord)
coordenadas<- as.data.frame(coordenadas)

# carregamento mcmc_funcion.R
mcmc_dados = syn_mcmc(dados_sem_coord, coordenadas, limits = c(0,10,0,10),
                      grid = 10, S = 1250, burn = 250,
                      continuous = 4,
                      spatial_beta = FALSE,
                      return_parameters = TRUE)

## Graficos

grid = 10
S = mcmc_dados[["S"]]
burn = mcmc_dados[["burn"]]
lambda = mcmc_dados[["lambda"]]
media.lambda = mcmc_dados[["media.lambda"]]
alfa = mcmc_dados[["alfa"]]
mu = mcmc_dados[["mu"]]
theta = mcmc_dados[["theta"]]
# tau.theta = mcmc_dados[["tau.theta"]]
phi = mcmc_dados[["phi"]]
# tau.phi = mcmc_dados[["tau.phi"]]
epsilon = mcmc_dados[["epsilon"]]
# tau.e = mcmc_dados[["tau.e"]]


#Gráficos de convergência

#traceplots dos parâmetros do modelo


plot_traceplots <- function(S, alfa, mu,thau.phi, tau.theta, tau.e ) {

  num_plots <- 0  # Variável para contar o número de gráficos
  plots <- list()  # Lista para armazenar os gráficos

  if (!is.null(mu)) {
    num_plots <- num_plots + 1
    plots[[num_plots]]<-plot(mu[1:S], type = "l", xlab = "Iterações", ylab = "mu")
  }

  if (!is.null(alfa)) {
    num_alfa <- nrow(alfa)
    for (i in 1:num_alfa) {
      num_plots <- num_plots + 1
      plots[[num_plots]]<-plot(alfa[i, 1:S], type = "l", xlab = "Iterações", ylab = paste("alfa (", i, ")"))
    }
  }

  if (!is.null(theta)) {
    num_theta <- nrow(theta)
    for (i in 1:num_theta) {
      num_plots <- num_plots + 1
      plots[[num_plots]]<-plot(theta[i, 1:S], type = "l", xlab = "Iterações", ylab = paste("theta (", i, ")")
    }
  }

  if (!is.null(phi)) {
    num_phi <- dim(phi)[]#?
    for (i in 1:num_phi) {
      num_plots <- num_plots + 1
      plots[[num_plots]]<-plots[[num_plots]]<-plot(phi[1:S, i], type = "l", xlab = "Iterações", ylab = paste("phi (", i, ")")
    }
  }

  if (!is.null(epsilon)) {
    num_epsilon <- dim(epsilon)[]#?
    for (i in 1:num_epsilon) {
      plots[[num_plots]]<-plot(epsilon[1:S, i], type = "l", xlab = "Iterações", ylab = paste("epsilon (", i, ")")
    }
  }

  # num_cols <- 3 #numero de colunas
  # num_rows <- ceiling(num_plots / num_cols) #numeros de linhas necessárias
  #
  # par(mfrow = c(num_rows, num_cols))

  # Plote os gráficos
  for (i in 1:num_plots) {
    plot(plots[[i]])
  }

  return(plots)  # Retorna a lista de gráficos
}



#plot de lambda
plot_lambda <- function(lambda, n_lambda, S_inicio, S) {
  num_combinacoes <- nrow(lambda)
  plot_graphs <- list()

  for (i in 1:num_combinacoes) {
    plot_graphs[[i]] <- plot(lambda[n_lambda, S_inicio:S, i], type = "l", xlab = "Iterações", ylab = paste("Lambda (", n_lambda, ") combinação ", i))
  }

  return(plot_graphs)
}

# Plot para avaliação do modelo
# para plotar gráficos de intensidade para diferentes combinações
plot_media_lambda <- function(media_lambda, grid, dados, comb) {
  corDegrade <- colorRampPalette(c("light yellow","red"))
  num_combinacoes <- length(media_lambda) #não sei qual o formato do media_lambda
  num_rows <- ceiling(num_combinacoes / 3)  # 3 gráficos por linha- celling arredonda para cima
  num_cols <- 3  # 3 gráficos por coluna

  par(mfrow = c(num_rows, num_cols))


  #nao ficou claro de onde tirou o lonvec e latvec (na duvida, vou classsificar como vetor)
  lonvec = dados$lon
  lacvec = dados$lat

  for (i in 1:num_combinacoes) {
    image(lonvec, latvec, matrix(media_lambda[, i], nrow = grid, ncol = grid), col = corDegrade(10),
          main = paste("Intensidade", colnames(media_lambda)[i]))
    points(dados$lon[comb == i], dados$lat[comb == i])
    rect(4, 4, 6, 5, border = "black")
  }

  dev.off()
}

