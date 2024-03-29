#' @title  Application of the Markov Chain Monte Carlo using the nimble package
#'
#' @name  syn_mcmc_nimble
#'
#' @description  To obtain the synthetic coordinates we are going to use the package nimble to apply the MCMC.
#' Function "syn_mcmc_nimble" receives our database uses variables that we got in the function \code{link:prepare_data} to do the mcmc with nimble.
#' By the end of this function we get the parameter \code{lambda} that will be required when creating the synthetic coordinates.
#'
#' @param  dataset   A data frame with all the information except the coordinates
#' @param  coord   An object with two columns indicating the latitude and longitude respectively of the elements in the dataset
#' @param  limits An object that is a vector of the dimensions where will be create the grids passed through the sequence of xmin, xmax, ymin, ymax. The default is create by using the maximum and the minimum of the coords object.
#' @param  grid  The grid represents the quantities of divisions that will be made in the location. Bigger the grid, closer the synthetic coordinates are to the real coordinates. With a default result of (grid = 10)
#' @param  nchain   Quantities of chains that will be made. With a default result of (nchain = 1)
#' @param  S   Quantities of simulations that will be made. With a default result of (S = 5000)
#' @param  burn   The number of simulations that will be burned to warm-up the \code{mcmc}. With a default result of (\code{burn} = 1000)
#' @param  continuous  Option so the user can warn the function for the presence of continuous variables in the dataset. You must indicate the columns numbers that contain continuous variables in a vector
#' @param  spatial_beta  Option so the user can choose to use a spatial \code{beta} parameter. You must input a vector containing which \code{beta} you want to add the spatiality. If you want to select all beta, you can use the logical argument TRUE
#' @param  return_parameters  Option to return the result of the parameters. With default not returning the parameters
#'
#' @return  Depending on the \code{return_parameters} parameter, this function can return only the \code{lambda} parameter or all other significant parameters too.
#'
#' @references
#' NUNES, Letícia. Métodos de Simulação de Dados Geográficos Sintéticos Para Bases Confidenciais. *Dissertação de Mestrado*, [s. l.], 2018.
#' Disponível em: \url:{http//est.ufmg.br/portal/arquivos/mestrado/dissertacoes/dissertacao_Leticia_Silva_Nunes.pdf}. Acesso em: 2 mar. 2022.
#'
#' @examples
#'   syn_mcmc_nimble(dataset = my_database, S = 2500, burn = 500, return_parameters = TRUE)
#'
#' @import nimble
#'
#' @export

syn_mcmc_nimble <- function(dataset, coord, limits = c(), grid = 10,
                     #nchains = 1, maybe edit later the code to suport multiple chains
                     S = 5000, burn = 1000,
                     continuous = FALSE, spatial_beta = FALSE,
                     return_parameters = FALSE){

  if (length(continuous) == 1)
    if (continuous == TRUE)
      stop("The continuous object must indicates which columns in the dataset correspond to continuous variables")


  if(is.numeric(continuous)){
    saida = prepare_data(dataset, coord, limits, grid, continuous)
  }else{
    saida = prepare_data(dataset, coord, limits, grid)
  }

  # Assigning the elements of the output to new objects

  mapply(assign, names(saida), saida, MoreArgs=list(envir = globalenv()))

  # Fitting additional variables for spatial modeling

  #Adjusting variables to recreate the matrix of combinations####
  library(Matrix)
  matrix.ind.a <- stack(setNames(ind.a, seq_along(ind.a)))
  matrix.ind.a <- as.matrix(sparseMatrix(as.numeric(matrix.ind.a[,2]), matrix.ind.a[,1], x=1))

  #Neighborhood function - maybe will be possible to specify your own in the future####
  neigh = cell2nb((grid),(grid), type="queen")
  nb.nimble <- nb2WB(neigh)

  #b_epsilon - an hyperparameter that soon make it possible for the user to change it####
  b_epsilon = 0.1

  #z.pad - variable necessary for modeling with a continuous variable
  if(exists('z.pad')){
    C <- length(z.pad[1,1,])


    z.pad[which(is.na(z.pad))]<-0

    if(C > 1){
      z <- matrix(NA,G*C,B)
      for(i in 1:G){
        for(j in 1:B){
          for(c in 1:C){
            z[i+G*(c-1),j] <- z.pad[i,j,c]
          }
        }
      }
    }else{
      z <- z.pad[,,1]
    }
  }else{
    C <- FALSE
  }

  # NimbleCode
  if(C == FALSE){
    #No continuous variable -> C = FALSE
    Nimble <- discrete_nimble(G, B, nx, matrix.ind.a, nb.nimble, ni, b_epsilon)
  }else if(C == 1){
    #Only one continuous variable -> C == 1
    if(spatial_beta == FALSE){
      Nimble <- single_nimble_beta(G, B, nx, matrix.ind.a, nb.nimble, ni, z, b_epsilon)
    }else{
      Nimble <- single_nimble_spatial_beta(G, B, nx, matrix.ind.a, nb.nimble, ni, z, b_epsilon)
    }
  }else{
    #Multiple continuous variaables -> C > 1
    if(spatial_beta == FALSE){
      Nimble <- mult_nimble_beta(G, B, nx, matrix.ind.a, nb.nimble, ni, z, C, b_epsilon)
    }else{
      Nimble <- mult_nimble_spatial_beta(G, B, nx, matrix.ind.a, nb.nimble, ni, z, C, b_epsilon)
    }
  }

  # Application of the MCMC
  if(return_parameters == FALSE){
    mcmc.out <- nimbleMCMC(code = Nimble$Code, constants = Nimble$Consts,
                         data = Nimble$Data, inits = Nimble$Inits,
                         nchains = 1, niter = S, nburnin = burn,
                         summary = TRUE, WAIC = TRUE,
                         monitors = c('lambda')
                         )
  }else if(return_parameters == TRUE & C == FALSE){
    mcmc.out <- nimbleMCMC(code = Nimble$Code, constants = Nimble$Consts,
                           data = Nimble$Data, inits = Nimble$Inits,
                           nchains = 1, niter = S, nburnin = burn,
                           summary = TRUE, WAIC = TRUE,
                           monitors = c('mu',
                                        'alfa',
                                        'theta',
                                        'tau_theta',
                                        'phi',
                                        'tau_phi',
                                        'epsilon',
                                        'tau_e',
                                        'lambda'
                           ))
  }else if(return_parameters == TRUE & C != FALSE){
    mcmc.out <- nimbleMCMC(code = Nimble$Code, constants = Nimble$Consts,
                           data = Nimble$Data, inits = Nimble$Inits,
                           nchains = 1, niter = S, nburnin = burn,
                           summary = TRUE, WAIC = TRUE,
                           monitors = c('mu',
                                        'alfa',
                                        'theta',
                                        'tau_theta',
                                        'phi',
                                        'tau_phi',
                                        'epsilon',
                                        'tau_e',
                                        'beta',
                                        'tau_beta',
                                        'lambda'
                           ))
  }

  # Organizing parameters to return at the end of the function
  sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,6)=="lambda")]; lambda <- array(,c(G,S-burn,B))
  for(j in 1:(S-burn)){
    for(k in 1:B){
      for(i in 1:G){lambda[i,j,k] <- sup[j,i+G*(k-1)]}
    }
  }

  media.lambda <- array(NA, c(G, B))
  for(i in 1:G){
    for(j in 1:B){media.lambda[i,j] <- mean(mcmc.out$samples[ , paste0("lambda[", i, ", ", j,"]")])}
  }

  if(return_parameters==TRUE){
    #alfa
    sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,4)=="alfa")]; alfa <- array(,c(sum(nx-1),S-burn))
    for(i in 1:(sum(nx-1))){
      for(j in 1:(S-burn)){alfa[i,j] <- sup[j,i]}
    }

    #mu
    sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,2)=="mu")]; mu <- array(,c(S-burn))
    for(i in 1:(S-burn)){mu[i] <- sup[i]}

    #theta
    sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,5)=="theta")]; theta <- array(,c(G,S-burn))
    for(i in 1:G){
      for(j in 1:(S-burn)){theta[i,j] <- sup[j,i]}
    }

    #tau.theta
    sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,9)=="tau_theta")]; tau.theta <- array(,c(S-burn))
    for(i in 1:(S-burn)){tau.theta[i] <- sup[i]}

    #phi
    sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,3)=="phi")]; phi <- array(,c(G,S-burn,sum(nx-1)))
    for(k in 1:(sum(nx-1))){
      for(j in 1:(S-burn)){
        for(i in 1:G){phi[i,j,k] <- sup[j,i+G*(k-1)]}
      }
    }

    #tau.phi
    sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,7)=="tau_phi")]; tau.phi <- array(,c(sum(nx-1),S-burn))
    for(i in 1:(sum(nx-1))){
      for(j in 1:(S-burn)){tau.phi[i,j] <- sup[j,i]}
    }

    #epsilon
    sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,7)=="epsilon")]; epsilon <- array(,c(G,S-burn,B))
    for(k in 1:B){
      for(j in 1:(S-burn)){
        for(i in 1:G){epsilon[i,j,k] <- sup[j,i+G*(k-1)]}
      }
    }

    #tau.e
    sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,5)=="tau_e")]; tau.e <- array(,c(S-burn))
    for(i in 1:(S-burn)){tau.e[i] <- sup[i]}

    if(C != FALSE){
      #beta
      sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,4)=="beta")]; beta <- array(,c(G,S-burn,C))
      for(k in 1:C){
        for(j in 1:(S-burn)){
          for(i in 1:G){beta[i,j,k] <- sup[j,i+G*(k-1)] }
        }
      }

      #tau.beta
      sup <- mcmc.out$samples[,which(substr(colnames(mcmc.out$samples),1,8)=="tau_beta")]; tau.beta <- array(,c(C, S-burn))
      for(i in 1:C){
        for(j in 1:(S-burn)){tau.beta[i,j] <- ifelse(C>1,sup[j,i],sup[j])}
      }
    }
  }

  if(return_parameters == FALSE){
    return(list(S=S, burn=burn, lambda=lambda, media.lambda=media.lambda))
  }else if(return_parameters == TRUE & C == FALSE){
    return(list(S=S, burn=burn, lambda=lambda, media.lambda=media.lambda,
                alfa=alfa, mu=mu, theta=theta, tau.theta=tau.theta, phi=phi, tau.phi=tau.phi, epsilon=epsilon, tau.e=tau.e))
  }else if(return_parameters == TRUE & C != FALSE){
    return(list(S=S, burn=burn, lambda=lambda, media.lambda=media.lambda,
                alfa=alfa, mu=mu, theta=theta, tau.theta=tau.theta, phi=phi, tau.phi=tau.phi, beta=beta, tau.beta=tau.beta, epsilon=epsilon, tau.e=tau.e))
  }

}
