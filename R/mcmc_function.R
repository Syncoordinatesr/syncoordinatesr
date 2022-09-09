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
#' @return  Depending on the \code{return_parameters} parameter, this function can return only the \code{lambda} parameter or all other significant parameters too.
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

  if(continuous != FALSE){
    saida = prepare_data(dataset, coord, grid, continuous)
  } else{
    saida = prepare_data(dataset, coord, grid)
  }

  # Assigning the elements of the output to new objects

  mapply(assign, names(saida), saida, MoreArgs=list(envir = globalenv()))

  acomb = function(x,i) i %in% x

  # Initializing the parameters

  alfa = matrix(0, sum(nx - 1), S)
  alfa[,1] = 0

  mu = numeric(S)
  mu[1] = 1

  # Adding beta
  if(continuous != FALSE){
    beta = array(NA, c(G, S, vZ))
    beta[,1,] = 1
    beta.atual = beta[,1,]
  }

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

  # Hyperparameters (soon make it possible for the user to change it)

  atheta = abeta = 0.1
  btheta = bbeta = 0.1
  aphi = 0.1
  bphi = 0.1
  vmu = 5
  valfa = 5
  vbeta = 5
  m.bar = mean(ni)
  be = 0.1; ae = m.bar*(0.7^2)*be
  #eta = array(NA, c(G, B, vZ))
  #eta.atual = array(NA, c(G, B, vZ))

  # Linear predictor

  #if(continuous != FALSE){
  #  for(k in 1:lenght(continuous)){
  #    for (j in 1:B){
  #      for (i in 1:G){
  #        eta[i,j,k] = mu[1] + sum(alfa[Z[,j],1]) + theta[i] + beta[i,1,k]%*%z.pad[i,j,k] + sum(phi[i,1,Z[,j]]) + epsilon[i,j]
  #        #eta.atual[i,j,k] = eta[i,j,k]
  #      }
  #    }
  #  }
  #}

  # ANTERIORMENTE NO DELA
  #for (j in 1:B){
  #  for (i in 1:G){
  #    eta[i,j]=mu[1]+sum(alfa[Z[,j],1])+theta[i]+beta[i]%*%z.pad[i,j]+sum(phi[i,1,Z[,j]])+epsilon[i,j]
  #    eta.atual[i,j]=eta[i,j]
  #  }
  #}

  #gama = array(data=0, dim=c(G, ifelse(continuous != FALSE, vZ, 1), B))
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
      if(continuous != FALSE){
        gama[,1,i] = gama[,1,i] + ifelse(vZ > 1,
                                         apply(beta[,1,]*z.pad[,i,], MAR=1, FUN=sum, na.rm = TRUE),
                                         ifelse(is.na(beta[,1,]*z.pad[,i,]), 0, beta[,1,]*z.pad[,i,]))
      }
    }
  gama.atual = gama[,1,]

  lambda = array(data=0,dim=c(G, S, B))
  lambda[,1,] = exp(gama)

  media.lambda = matrix(0, G, B)

  ##########################
  # Start of the mcmc

  controle=MfU.Control(n=1,slice.w=.01,slice.m=10000,slice.lower=0,slice.upper=1)

  for (k in 2:S){ #Traduzir pro nimble

    cat("MCMC simulation ",k," of ",S,"\n")

    gama = gama.atual - mu[k-1]
    sumeta = sum(exp(gama))

    # Estimation of mu

    mu[k] = ars(1, muf, mufprima,
                lb=T, xlb=-100, ub=T, xub=100,
                sumeta=sumeta, vmu=vmu, ci_b=n)
    gama.atual = gama + mu[k]

    # Estimation of alfa

    if(continuous != FALSE){

      for(t in 1:dim(alfa)[1]){
        temp = apply(Z,MAR=2,FUN=acomb,t) #Checar se realmente o interesse é na transposta de Z
        n.alfa=sum(ci_b[,temp])
        gama=gama.atual[,temp]-alfa[t,(k-1)]
        sumeta.a=sum(exp(gama))
        alfa[t,k] <- ars(1, muf,mufprima,lb=T,xlb=-100,ub=T,xub=100,sumeta=sumeta.a,vmu=valfa,ci_b=n.alfa)
        gama.atual[,temp] = gama + alfa[t,k]
      }
    } else{
      for(t in 1:dim(alfa)[1]){
        # n.alfa = sum(ci_b[sub.a[[t]],]) # original code by Leticia
        ## it's accessing the wrong dimension of ci_b, it should be:
        n.alfa = sum(ci_b[ , sub.a[[t]]]) # changed by Thais - Feb/22
        gama = gama.atual[,sub.a[[t]]] - alfa[t,(k-1)]
        sumeta.a = sum(exp(gama))
        alfa[t,k] <- ars(1, muf, mufprima,
                         lb=T, xlb=-100, ub=T, xub=100,
                         sumeta=sumeta.a, vmu=valfa, ci_b=n.alfa)
        gama.atual[,sub.a[[t]]] = gama + alfa[t,k]
      }
    }

    # Estimation of theta

    # n.theta = apply(ci_b, MAR=2, FUN=sum) ## original code by Leticia
    ## it was returning object with wrong dimension, and causing errors in the 'ars' call bellow - I am assuming this should be sum(ci_b) by b=1,...,B
    n.theta = apply(ci_b, MAR=1, FUN=sum) ## changed by Thais (Feb/2022)
    gama = gama.atual - theta[,(k-1)]
    for(g in 1:G){
      sum.theta = sum(exp(gama[g,]))
      bar = W[g,]%*%theta.atual/ni[g]
      theta.atual[g] = ars(1, thetaf, thetafprima, ns=1000,
                           lb=T, xlb=-100, ub=T, xub=100, ci_b=n.theta[g],
                           sumeta=sum.theta, ni=ni[g],
                           tau.f=tau.theta[k-1], bar.f=bar)
    }
    theta[,k] = theta.atual - sum(theta.atual)/G
    gama.atual = gama + theta[,k]

    # Estimation of phi - ATTENTION TO Z - Checar o que é o temp

    if(continuous != FALSE){

      for(t in 1:dim(phi)[3]){
        temp = apply(Z,MAR=2,FUN=acomb,t)
        n.phi = apply(ci_b[,temp],MAR=1,FUN=sum)
        gama = gama.atual[,temp] - phi[,(k-1),t]
        for(g in 1:G){
          sum.phi = sum(exp(gama[g,]))
          bar = W[g,]%*%phi.atual[,t]/ni[g]
          phi.atual[g,t] = ars(1,thetaf,thetafprima,lb=T,ns=1000,xlb=-100,ub=T,xub=100,ci_b=n.phi[g],sumeta=sum.phi,ni=ni[g],tau.f=tau.phi[t,(k-1)],bar.f=bar)
        }
        ## Sum equal to zero
        phi[,k,t] = phi.atual[,t] - sum(phi.atual[,t])/G
        gama.atual[,temp] = gama + phi[,k,t]
      }
    } else{

      for(t in 1:dim(phi)[3]){
        # n.phi = apply(saida$ci_b[sub.a[[t]],],MAR=2,FUN=sum) ## original code by Leticia, the output dimension seems wrong
        ## there is also a problem with the dimensions of ci_b
        n.phi = apply(ci_b[,sub.a[[t]]],MAR=1,FUN=sum) ## changed by Thais (Feb/2022)
        gama = gama.atual[,sub.a[[t]]] - phi[,(k-1),t]
        for(g in 1:G){
          sum.phi = sum(exp(gama[g,]))
          bar = W[g,]%*%phi.atual[,t]/ni[g]
          phi.atual[g,t] = ars(1, thetaf, thetafprima, ns=1000,
                               lb=T, xlb=-100, ub=T, xub=100,
                               ci_b=n.phi[g], sumeta=sum.phi, ni=ni[g],
                               tau.f=tau.phi[t,(k-1)], bar.f=bar)
        }
        phi[,k,t] = phi.atual[,t] - sum(phi.atual[,t])/G
        # gama.atual[,sub.a[[k]]] = eta + phi[,k,t] ## original by Leticia
        ## seems like it was saving the updated gama in the wrong positions
        ## k is the simulations' index, should not be used to save gama!
        gama.atual[,sub.a[[t]]] = gama + phi[,k,t] ## changed by Thais (Feb/22)
      }
    }

    ################################################################################################################################################
    #Adding beta - check it
    if(continuous != FALSE){
      #gama = gama.atual - beta[,(k-1)]*z.pad   <----  ANTES NO DELA
      for(i in 1:vZ){

        gama = gama.atual - ifelse(vZ > 1,
                                    apply(beta[,(k-1),]*z.pad[,1,], MAR=1, FUN=sum),
                                     beta[,(k-1),]*z.pad[,1,])

        c.f = array(NA, c(G, B, vZ))
        c.f[,,i] = apply(ci_b*z.pad[,,i],MAR=1,FUN=sum) #CHECAR DIMENSÃO DO Z.PAD - FAZER UM LOOP PARA QUANTOS Z TIVER - cib sempre o mesmo inedependente do z
        for(g in 1:G){
          zib.vec = z.pad[g,,i] #CHECAR DIMENSÃO DO Z.PAD
          u <- 0
          u <- (MfU.Sample(x=logit(beta[g,k-1,i]),f=betaf,uni.sampler="slice",c.f=c.f[g+((i-1)*(G*B))],gama=gama[g,],vbeta=vbeta,zib.vec=zib.vec,control=controle)) #ANTES c.f=c.f[g]
          beta[g,k,i] = inv.logit(u) #CHECAR DIMENSÃO DO BETA - aumentar
        }
        #if(spatial_beta == FALSE){
        #  gama.atual = gama + beta[,k]*z.pad
        #} else{
        #  beta[,k] = beta.atual - sum(beta.atual)/G
        #  gama.atual = gama + beta[,k]*z.pad
        #}
      }
    }

    # Estimation of epsilon

    for(g in 1:G){
      for(j in 1:B){
        gama = gama.atual[g,j] - epsilon[g,1,j]
        sum.aux = exp(gama)
        epsilon[g,1,j] = ars(1, ef, efprima,
                             lb=T, xlb=-100, ub=T, xub=100,
                             #ci_b=ci_b[j,g], sumeta=sum.aux, ## original code by Leticia - dimensions were exchanges
                             ci_b=ci_b[g,j], sumeta=sum.aux, ## code changed by Thais (Feb/22)
                             tau.e=tau.e[k-1])
        gama.atual[g,j] = gama + epsilon[g,1,j]
      }
    }

    # compute the up-to-date mean of lambda, in case of returning just the mean
    if(k > burn){
      media.lambda = (media.lambda*((k-burn-1)/(k-burn))) +
        # (exp(eta.atual)/(k-burn)) # original code by Leticia ??
        (exp(gama.atual)/(k-burn)) # changed by Thais (Feb/2022)
    }

    sum.theta = t(theta[,k])%*%(diag(ni)-W)%*%(theta[,k])
    tau.theta[k] = rgamma(1,atheta+n/2,rate=btheta+sum.theta/2)

    for(m in 1:dim(tau.phi)[1]){
      sum.phi = t(phi[,k,m])%*%(diag(ni)-W)%*%(phi[,k,m])
      tau.phi[m,k] = rgamma(1,aphi+G/2,rate=bphi+sum.phi/2)
    }

    ## Update of tau.e
    sum.e = sum(epsilon^2)
    tau.e[k] = rgamma(1,ae+(G*B)/2,rate=be+sum.e/2)

    if(continuous != FALSE){
      sum.beta = t(beta[,k,]%*%(diag(ni)-W)%*%beta[,k,])
      tau.beta[k] = rgamma(1,abeta+n/2,rate=bbeta+sum.beta/2)
    }

    # what if we want to return the entire lambda?? CHECK SIZE! - Thais
    lambda[,k,] = exp(gama.atual)

  } #mcmc end

  if (return_parameters == FALSE){
    return(list(S=S, burn=burn, lambda=lambda, media.lambda=media.lambda))
  }
  else{
    return(list(S=S, burn=burn, lambda=lambda, media.lambda=media.lambda,
                alfa=alfa, mu=mu, theta=theta, tau.theta=tau.theta, phi=phi, tau.phi=tau.phi, epsilon=epsilon, tau.e=tau.e))
  }
}
