#' @title Markov Chain Monte Carlo function
#'
#' @name  syn_mcmc
#'
#' @description  A MCMC function that will be useful to obtain the synthetic coordinates.
#' Function "syn_mcmc" receives our database uses variables that we got in the function \code{link:prepare_data} to do the mcmc.
#' By the end of this function we get the parameter \code{lambda} that will be required when creating the synthetic coordinates.
#'
#' @param  dataset   A data frame with all the information except the coordinates
#' @param  S   Quantities of simulations that will be made. With a default result of (S = 5000)
#' @param  burn   The number of simulations that will be burned to warm-up the mcmc. With a default result of (burn = 1000)
#' @param  continuous  Option so the user can warn the function for the presence of continuous variables in the dataset. With default considering only discrete variables
#' @param  spatial_beta  Option so you can choose to use a spatial beta parameter
#' @param  return_paramenters  Option to return the result of the parameters. With default not returning the parameters
#'
#' @return  Depending on the \code{return_parameters} parameter, this function can return only the \code{lambda} parameter, or all other significant parameters too.
#'
#' @references
#' NUNES, Letícia. Métodos de Simulação de Dados Geográficos Sintéticos Para Bases Confidenciais. *Dissertação de Mestrado*, [s. l.], 2018.
#' Disponível em: \url:{http//est.ufmg.br/portal/arquivos/mestrado/dissertacoes/dissertacao_Leticia_Silva_Nunes.pdf}. Acesso em: 2 mar. 2022.
#'
#' @examples
#'   syn_mcmc(dataset = my_database, S = 2500, burn = 500, return_parameters = TRUE)
#'
#' @import ars MfUSampler
#'
#' @export

syn_mcmc <- function(dataset, coord, grid = 10,
                          S = 5000, burn = 1000,
                          continuous = FALSE, spatial_beta = FALSE,
                          return_parameters = FALSE){

  if(continuous = TRUE){
    saida = prepare_data(dataset, coord, grid, continuous = TRUE)
  } else{
    saida = prepare_data(dataset, coord, grid)
  }

  # assigning the elements of the output to new objects
  mapply(assign, names(saida), saida, MoreArgs=list(envir = globalenv()))

  #####Não está encontrando as funções que estão na pasta Utils

  # Falta Z e z.pad

  acomb = function(x,i) i %in% x #Added

  # Inicializando os parâmetros

  alfa = matrix(0, sum(nx - 1), S)
  alfa[,1] = 0

  mu = numeric(S)
  mu[1] = 1

  # Adding beta
  beta = matrix(0,G,S)
  beta[,1] = 0
  beta.atual=beta[,1] #Added

  theta = matrix(0, G, S)
  theta[,1] = theta[,1] - sum(theta[,1])/(G)
  theta.atual = theta[,1]

  tau.theta = numeric(S)
  tau.theta[1] = 1

  phi = array(data=0,dim=c(G,S,sum(nx-1)))
  for(i in 1:dim(phi)[3]){
    phi[,1,i] = phi[,1,i] - mean(phi[,1,i])
  }
  phi.atual = phi[,1,]

  tau.phi = matrix(1,dim(phi)[3],S)
  tau.phi[,1] = 1

  tau.beta = numeric(S)
  tau.beta[1] = 1

  epsilon = array(data=0, dim=c(G, 1, B))

  tau.e = numeric(S)
  tau.e[1] = 1

  # Hyperparameters (no futuro tornar possível para o usuário alterar)

  atheta = abeta = 0.1 #Changed
  btheta = bbeta = 0.1 #Changed
  aphi = 0.1
  bphi = 0.1
  vmu = 5
  valfa = 5
  vbeta = 5 #Added
  m.bar = mean(ni)
  be = 0.1; ae = m.bar*(0.7^2)*be
  eta = matrix(NA,G,B) #Added
  eta.atual = matrix(NA,G,B) #Added

  # Preditor linear

  gama = array(data=0, dim=c(G, 1, B))
  for(i in b){
    if(length(ind.a[[i]])==0){
      gama[,1,i] = log(n) + mu[1] + theta[,1] + epsilon[,1,i]
    }
    else{
      gama[,1,i] = log(n) + mu[1] + sum(alfa[ind.a[[i]],1]) +
        theta[,1] + ifelse(length(ind.a[[i]])>1,
                           apply(phi[,1,ind.a[[i]]],MAR=1,FUN=sum),
                           phi[,1,ind.a[[i]]]) + epsilon[,1,i]
    }
  }
  gama.atual = gama[,1,]

  lambda = array(data=0,dim=c(G, S, B))

  media.lambda = matrix(0, G, B)

  ##########################
  # começo real do mcmc

  controle=MfU.Control(n=1,slice.w=.01,slice.m=10000,slice.lower=0,slice.upper=1)

  for (k in 2:S){

    cat("MCMC simulation ",k," of ",S,"\n")

    eta = gama.atual - mu[k-1]
    sumeta = sum(exp(eta))
    mu[k] = ars(1, muf, mufprima,
                lb=T, xlb=-100, ub=T, xub=100,
                sumeta=sumeta, vmu=vmu, ci_b=saida$n)
    gama.atual = eta + mu[k]

    if(continuous == FALSE){

      for(t in 1:dim(alfa)[1]){
        # n.alfa = sum(ci_b[sub.a[[t]],]) # original code by Leticia
        ## it's accessing the wrong dimension of ci_b, it should be:
        n.alfa = sum(ci_b[ , sub.a[[t]]]) # changed by Thais - Feb/22
        eta = gama.atual[,sub.a[[t]]] - alfa[t,(k-1)]
        sumeta.a = sum(exp(eta))
        alfa[t,k] <- ars(1, muf, mufprima,
                         lb=T, xlb=-100, ub=T, xub=100,
                         sumeta=sumeta.a, vmu=valfa, ci_b=n.alfa)
        gama.atual[,sub.a[[t]]] = eta + alfa[t,k]
      }
    } else{
      for(t in 1:dim(alfa)[1]){
        temp = apply(Z,MAR=2,FUN=acomb,t)
        n.alfa=sum(ci_b[,temp])
        eta=gama.atual[,temp]-alfa[t,(k-1)] #Usando gama.atual ao inves de eta.atual
        sumeta.a=sum(exp(eta))
        alfa[t,k] <- ars(1, muf,mufprima,lb=T,xlb=-100,ub=T,xub=100,sumeta=sumeta.a,vmu=valfa,ci_b=n.alfa)
        gama.atual[,temp] = eta + alfa[t,k] #Usando gama.atual ao inves de eta.atual
      }
    }

    # n.theta = apply(ci_b, MAR=2, FUN=sum) ## original code by Leticia
    ## it was returning object with wrong dimension, and causing errors in the 'ars' call bellow - I am assuming this should be sum(ci_b) by b=1,...,B
    n.theta = apply(ci_b, MAR=1, FUN=sum) ## changed by Thais (Feb/2022)
    eta = gama.atual - theta[,(k-1)]
    for(g in 1:saida$G){
      sum.theta = sum(exp(eta[g,]))
      bar = W[g,]%*%theta.atual/ni[g]
      theta.atual[g] = ars(1, thetaf, thetafprima, ns=1000,
                           lb=T, xlb=-100, ub=T, xub=100, ci_b=n.theta[g],
                           sumeta=sum.theta, ni=ni[g],
                           tau.f=tau.theta[k-1], bar.f=bar)
    }
    theta[,k] = theta.atual - sum(theta.atual)/saida$G
    gama.atual = eta + theta[,k]

    if(continuous == FALSE){

      for(t in 1:dim(phi)[3]){
        # n.phi = apply(saida$ci_b[sub.a[[t]],],MAR=2,FUN=sum) ## original code by Leticia, the output dimension seems wrong
        ## there is also a problem with the dimensions of ci_b
        n.phi = apply(saida$ci_b[,sub.a[[t]]],MAR=1,FUN=sum) ## changed by Thais (Feb/2022)
        eta = gama.atual[,sub.a[[t]]] - phi[,(k-1),t]
        for(g in 1:saida$G){
          sum.phi = sum(exp(eta[g,]))
          bar = W[g,]%*%phi.atual[,t]/ni[g]
          phi.atual[g,t] = ars(1, thetaf, thetafprima, ns=1000,
                               lb=T, xlb=-100, ub=T, xub=100,
                               ci_b=n.phi[g], sumeta=sum.phi, ni=ni[g],
                               tau.f=tau.phi[t,(k-1)], bar.f=bar)
        }
        phi[,k,t] = phi.atual[,t] - sum(phi.atual[,t])/saida$G
        # gama.atual[,sub.a[[k]]] = eta + phi[,k,t] ## original by Leticia
        ## seems like it was saving the updated gama in the wrong positions
        ## k is the simulations' index, should not be used to save gama!
        gama.atual[,sub.a[[t]]] = eta + phi[,k,t] ## changed by Thais (Feb/22)
      }
    } else{
      for(t in 1:dim(phi)[3]){
        temp = apply(Z,MAR=2,FUN=acomb,t)
        n.phi = apply(ci_b[,temp],MAR=1,FUN=sum)
        eta = gama.atual[,temp] - phi[,(k-1),t] #Usando gama.atual ao inves de eta.atual
        for(g in 1:G){
          sum.phi = sum(exp(eta[g,]))
          bar = W[g,]%*%phi.atual[,t]/ni[g]
          phi.atual[g,t] = ars(1,thetaf,thetafprima,lb=T,ns=1000,xlb=-100,ub=T,xub=100,ci_b=n.phi[g],sumeta=sum.phi,ni=ni[g],tau.f=tau.phi[t,(k-1)],bar.f=bar)
        }
        ## Soma igual a zero
        phi[,k,t] = phi.atual[,t] - sum(phi.atual[,t])/G
        gama.atual[,temp] = eta + phi[,k,t] #Usando gama.atual ao inves de eta.atual
      }
    }

    #Adding beta
    eta = gama.atual - beta[,(k-1)]*z.pad #Usando gama.atual ao inves de eta.atual ######ERRO
    c.f = apply(ci_b*z.pad,MAR=1,FUN=sum)
    for(g in 1:G){
      zib.vec = z.pad[g,]
      u <- 0
      u <- (MfU.Sample(x=logit(beta[g,k-1]),f=betaf,uni.sampler="slice",c.f=c.f[g],eta=eta[g,],vbeta=vbeta,zib.vec=zib.vec,control=controle))
      beta[g,k] = inv.logit(u)
    }
    if(spatial_beta == FALSE){
      gama.atual = eta + beta[,k]*z.pad #Usando gama.atual ao inves de eta.atual
    } else{
      beta[,k] = beta.atual - sum(beta.atual)/G
      gama.atual = eta + beta[,k]*z.pad #Usando gama.atual ao inves de eta.atual
    }

    for(g in 1:G){ #Changing saida$G for G
      for(j in 1:B){ #Changing saida$B for B
        eta = gama.atual[g,j] - epsilon[g,1,j]
        sum.aux = exp(eta)
        epsilon[g,1,j] = ars(1, ef, efprima,
                             lb=T, xlb=-100, ub=T, xub=100,
                             #ci_b=ci_b[j,g], sumeta=sum.aux, ## original code by Leticia - dimensions were exchanges
                             ci_b=ci_b[g,j], sumeta=sum.aux, ## code changed by Thais (Feb/22)
                             tau.e=tau.e[k-1])
        gama.atual[g,j] = eta + epsilon[g,1,j]
      }
    }

    # compute the up-to-date mean of lambda, in case of returning just the mean
    if(k > burn){
      media.lambda = (media.lambda*((k-burn-1)/(k-burn))) +
        # (exp(eta.atual)/(k-burn)) # original code by Leticia ??
        (exp(gama.atual)/(k-burn)) # changed by Thais (Feb/2022)
    }

    sum.theta = t(theta[,k])%*%(diag(ni)-W)%*%(theta[,k])
    tau.theta[k] = rgamma(1,atheta+saida$n/2,rate=btheta+sum.theta/2)

    for(m in 1:dim(tau.phi)[1]){
      sum.phi = t(phi[,k,m])%*%(diag(ni)-W)%*%(phi[,k,m])
      tau.phi[m,k] = rgamma(1,aphi+saida$G/2,rate=bphi+sum.phi/2)

      ## atualizando tau.e
      sum.e = sum(epsilon^2)
      tau.e[k] = rgamma(1,ae+(saida$G*saida$B)/2,rate=be+sum.e/2)
    }

    # what if we want to return the entire lambda?? CHECK SIZE! - Thais
    lambda[,k,] = exp(gama.atual)

  }

  if (return_parameters == FALSE){
    return(list(S=S, burn=burn, lambda=lambda, media.lambda=media.lambda))
  }
  else{
    return(list(S=S, burn=burn, lambda=lambda, media.lambda=media.lambda,
                alfa=alfa, mu=mu, theta=theta, tau.theta=tau.theta, phi=phi, tau.phi=tau.phi, epsilon=epsilon, tau.e=tau.e))
  }


}
