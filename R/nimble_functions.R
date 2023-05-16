#Caso discreto####
discrete_nimble <- function(G, B, nx, matrix.ind.a, nb.nimble, ni, b_epsilon){
  pumpCode <- nimbleCode({
    #Calculation of alfa_b for each combination
    for (j in 1:B){
      for(k in 1:a){
        sup_alfa[j,k] <- alfa[k] * matrix.ind.a[j,k]
      }

      alfa_b[j] <- sum(sup_alfa[j,1:a])
    }

    #Calculation of phi_b for each combination
    for(i in 1:G){
      for (j in 1:B){
        for(k in 1:a){
          sup_phi[i,j,k] <- phi[i,k] * matrix.ind.a[j,k]
        }

        phi_b[i,j] <- sum(sup_phi[i,j,1:a])
      }
    }

    #Modeling lambda and ci_b
    for(i in 1:G){
      for (j in 1:B){
        log(lambda[i,j]) <-  mu +
          alfa_b[j] +
          theta[i] +
          phi_b[i,j] +
          epsilon[i,j]# + n
        ci_b[i,j] ~ dpois(lambda[i,j])
      }
    }

    #Priori functions

    for(k in 1:a){
      phi[1:G,k] ~ dcar_normal(adj = w[1:L],
                               num = Nw[1:G],
                               tau = tau_phi[k],
                               zero_mean = 1)
      tau_phi[k] ~ dgamma(a_phi,b_phi)
    }


    theta[1:G] ~ dcar_normal(adj = w[1:L],
                             num = Nw[1:G],
                             tau = tau_theta,
                             zero_mean = 1)
    tau_theta ~ dgamma(a_theta,b_theta)


    for(k in 1:a){
      alfa[k] ~ dnorm(0, var = v_alfa_kj)
    }


    mu ~ dnorm(0, var = v_mu)


    for(i in 1:G){
      for(j in 1:B){
        epsilon[i,j] ~ dnorm(0, tau = tau_e)
      }
    }
    tau_e ~ dgamma(a_epsilon, b_epsilon)
  })

  pumpConsts <- list(
    #Model constants
    G = G, B = B,# n = n,                                          #generic
    a  = sum(nx-1), matrix.ind.a = matrix.ind.a,                  #combinations
    Nw = ni, w = nb.nimble$adj, L = sum(ni),                      #spatial
    #Priori functions constants
    v_mu = 5,                                                     #generic  constants
    v_alfa_kj = 5,                                                #combinations
    a_theta = 0.1, b_theta = 0.1,                                 #spatial constants
    a_phi = 0.1, b_phi = 0.1,                                     #spatial
    a_epsilon = mean(ni)*(0.7^2)*b_epsilon, b_epsilon = b_epsilon #error
  )
  pumpData <- list(ci_b = ci_b[1:pumpConsts$G, 1:pumpConsts$B])
  pumpInits <- list(
    #Initial values for the parameters
    mu = 1,                                                       #generic constants
    alfa = rep(0,pumpConsts$a),                                   #combinations
    theta = rep(1, pumpConsts$G),                                 #spatial constants
    phi = matrix(0, pumpConsts$G, pumpConsts$a),                  #spatial
    epsilon = matrix(0, pumpConsts$G, pumpConsts$B)               #error
  )

  return(list(Code = pumpCode,
              Consts = pumpConsts,
              Data = pumpData,
              Inits = pumpInits))
}

#Caso contínuo single normal beta####
single_nimble_beta <- function(G, B, nx, matrix.ind.a, nb.nimble, ni, z, b_epsilon){
  pumpCode <- nimbleCode({
    #Calculation of alfa_b for each combination####
    for (j in 1:B){
      for(k in 1:a){
        sup_alfa[j,k] <- alfa[k] * matrix.ind.a[j,k]
      }

      alfa_b[j] <- sum(sup_alfa[j,1:a])
    }

    #Calculation of phi_b for each combination####
    for(i in 1:G){
      for (j in 1:B){
        for(k in 1:a){
          sup_phi[i,j,k] <- phi[i,k] * matrix.ind.a[j,k]
        }

        phi_b[i,j] <- sum(sup_phi[i,j,1:a])
      }
    }

    #Modeling lambda and ci_b####
    for(i in 1:G){
      for (j in 1:B){
        log(lambda[i,j]) <-  mu +
          alfa_b[j] +
          theta[i] +
          phi_b[i,j] +
          beta[i]*z_pad[i,j] +
          epsilon[i,j]# + n
        ci_b[i,j] ~ dpois(lambda[i,j])
      }
    }

    #Priori functions####

    for(k in 1:a){
      phi[1:G,k] ~ dcar_normal(adj = w[1:L],
                               num = Nw[1:G],
                               tau = tau_phi[k],
                               zero_mean = 1)
      tau_phi[k] ~ dgamma(a_phi,b_phi)
    }


    theta[1:G] ~ dcar_normal(adj = w[1:L],
                             num = Nw[1:G],
                             tau = tau_theta,
                             zero_mean = 1)
    tau_theta ~ dgamma(a_theta,b_theta)


    for(k in 1:a){
      alfa[k] ~ dnorm(0, var = v_alfa_kj)
    }


    mu ~ dnorm(0, var = v_mu)


    for(i in 1:G){
      for(j in 1:B){
        epsilon[i,j] ~ dnorm(0, tau = tau_e)
      }
    }
    tau_e ~ dgamma(a_epsilon, b_epsilon)

    for(i in 1:G){
      beta[i] ~ dnorm(0, tau = tau_beta)
    }
    tau_beta ~ dgamma(a_beta, b_beta)

  })

  pumpConsts <- list(
    #Model constants
    G = G, B = B, #n = n,                                          #generic
    a  = sum(nx-1), matrix.ind.a = matrix.ind.a,                  #combinations
    Nw = ni, w = nb.nimble$adj, L = sum(ni),                      #spatial
    z_pad = z, #C = C,                       #continuous
    #Priori functions constants
    v_mu = 5,                                                     #generic  constants
    v_alfa_kj = 5,                                                #combinations
    a_theta = 0.1, b_theta = 0.1,                                 #spatial constants
    a_phi = 0.1, b_phi = 0.1,                                     #spatial
    a_epsilon = mean(ni)*(0.7^2)*b_epsilon, b_epsilon = b_epsilon, #error
    a_beta = 0.1, b_beta = 0.1                                    #continuous
  )
  pumpData <- list(ci_b = ci_b[1:pumpConsts$G, 1:pumpConsts$B])
  pumpInits <- list(
    #Initial values for the parameters
    mu = 1,                                                       #generic constants
    alfa = rep(0,pumpConsts$a),                                   #combinations
    theta = rep(1, pumpConsts$G),                                 #spatial constants
    phi = matrix(0, pumpConsts$G, pumpConsts$a),                  #spatial
    epsilon = matrix(0, pumpConsts$G, pumpConsts$B),              #error
    beta = rep(1, pumpConsts$G)                  #continuous
  )

  return(list(Code = pumpCode,
              Consts = pumpConsts,
              Data = pumpData,
              Inits = pumpInits))
}

#Caso contínuo single spatial beta####
single_nimble_beta_spatial <- function(G, B, nx, matrix.ind.a, nb.nimble, ni, z, b_epsilon){
  pumpCode <- nimbleCode({
    #Calculation of alfa_b for each combination
    for (j in 1:B){
      for(k in 1:a){
        sup_alfa[j,k] <- alfa[k] * matrix.ind.a[j,k]
      }

      alfa_b[j] <- sum(sup_alfa[j,1:a])
    }

    #Calculation of phi_b for each combination
    for(i in 1:G){
      for (j in 1:B){
        for(k in 1:a){
          sup_phi[i,j,k] <- phi[i,k] * matrix.ind.a[j,k]
        }

        phi_b[i,j] <- sum(sup_phi[i,j,1:a])
      }
    }

    #Modeling lambda and ci_b
    for(i in 1:G){
      for(j in 1:B){
        log(lambda[i,j]) <-  mu +
          alfa_b[j] +
          theta[i] +
          phi_b[i,j] +
          beta[i]*z_pad[i,j] +
          epsilon[i,j]# + n
        ci_b[i,j] ~ dpois(lambda[i,j])
      }
    }

    #Priori functions

    for(k in 1:a){
      phi[1:G,k] ~ dcar_normal(adj = w[1:L],
                               num = Nw[1:G],
                               tau = tau_phi[k],
                               zero_mean = 1)
      tau_phi[k] ~ dgamma(a_phi,b_phi)
    }


    theta[1:G] ~ dcar_normal(adj = w[1:L],
                             num = Nw[1:G],
                             tau = tau_theta,
                             zero_mean = 1)
    tau_theta ~ dgamma(a_theta,b_theta)


    for(k in 1:a){
      alfa[k] ~ dnorm(0, var = v_alfa_kj)
    }


    mu ~ dnorm(0, var = v_mu)


    for(i in 1:G){
      for(j in 1:B){
        epsilon[i,j] ~ dnorm(0, tau = tau_e)
      }
    }
    tau_e ~ dgamma(a_epsilon, b_epsilon)

    beta[1:G] ~ dcar_normal(adj = w[1:L],
                            num = Nw[1:G],
                            tau = tau_beta,
                            zero_mean = 1)
    tau_beta ~ dgamma(a_beta, b_beta)

  })

  pumpConsts <- list(
    #Model constants
    G = G, B = B,# n = n,                                         #generic
    a  = sum(nx-1), matrix.ind.a = matrix.ind.a,                  #combinations
    Nw = ni, w = nb.nimble$adj, L = sum(ni),                      #spatial
    z_pad = z, #C = C,                       #continuous
    #Priori functions constants
    v_mu = 5,                                                     #generic  constants
    v_alfa_kj = 5,                                                #combinations
    a_theta = 0.1, b_theta = 0.1,                                 #spatial constants
    a_phi = 0.1, b_phi = 0.1,                                     #spatial
    a_epsilon = mean(ni)*(0.7^2)*b_epsilon, b_epsilon = b_epsilon, #error
    a_beta = 0.1, b_beta = 0.1                                    #continuous
  )
  pumpData <- list(ci_b = ci_b[1:pumpConsts$G, 1:pumpConsts$B])
  pumpInits <- list(
    #Initial values for the parameters
    mu = 1,                                                       #generic constants
    alfa = rep(0,pumpConsts$a),                                   #combinations
    theta = rep(1, pumpConsts$G),                                 #spatial constants
    phi = matrix(0, pumpConsts$G, pumpConsts$a),                  #spatial
    epsilon = matrix(0, pumpConsts$G, pumpConsts$B),              #error
    beta = rep(1, pumpConsts$G)                 #continuous
  )

  return(list(Code = pumpCode,
              Consts = pumpConsts,
              Data = pumpData,
              Inits = pumpInits))
}

#Caso contínuo multiple normal beta####
mult_nimble_beta <- function(G, B, nx, matrix.ind.a, nb.nimble, ni, z, C, b_epsilon){
  pumpCode <- nimbleCode({
    #Calculation of alfa_b for each combination####
    for (j in 1:B){
      for(k in 1:a){
        sup_alfa[j,k] <- alfa[k] * matrix.ind.a[j,k]
      }

      alfa_b[j] <- sum(sup_alfa[j,1:a])
    }

    #Calculation of phi_b for each combination####
    for(i in 1:G){
      for (j in 1:B){
        for(k in 1:a){
          sup_phi[i,j,k] <- phi[i,k] * matrix.ind.a[j,k]
        }

        phi_b[i,j] <- sum(sup_phi[i,j,1:a])
      }
    }

    #Sum of the ponderated betas####
    for(i in 1:G){
      for (j in 1:B){
        for (c in 1:C){
          sup_beta[i,j,c] <- beta[i,c]*z_pad[(i+G*(c-1)),j]
        }
        beta_soma[i,j] <- sum(sup_beta[i,j,1:C])
      }
    }

    #Modeling lambda and ci_b####
    for(i in 1:G){
      for (j in 1:B){
        log(lambda[i,j]) <-  mu +
          alfa_b[j] +
          theta[i] +
          phi_b[i,j] +
          beta_soma[i,j] +
          epsilon[i,j]# + log(n)#se n=1 nao deveria mudar o resultado
        ci_b[i,j] ~ dpois(lambda[i,j])#e se multiplicar aqui por n
      }
    }

    #Priori functions####

    for(k in 1:a){
      phi[1:G,k] ~ dcar_normal(adj = w[1:L],
                               num = Nw[1:G],
                               tau = tau_phi[k],
                               zero_mean = 1)
      tau_phi[k] ~ dgamma(a_phi,b_phi)
    }


    theta[1:G] ~ dcar_normal(adj = w[1:L],
                             num = Nw[1:G],
                             tau = tau_theta,
                             zero_mean = 1)
    tau_theta ~ dgamma(a_theta,b_theta)


    for(k in 1:a){
      alfa[k] ~ dnorm(0, var = v_alfa_kj)
    }


    mu ~ dnorm(0, var = v_mu)


    for(i in 1:G){
      for(j in 1:B){
        epsilon[i,j] ~ dnorm(0, tau = tau_e)
      }
    }
    tau_e ~ dgamma(a_epsilon, b_epsilon)

    for(c in 1:C){
      tau_beta[c] ~ dgamma(a_beta, b_beta)

      for(i in 1:G){
        beta[i, c] ~ dnorm(0, tau = tau_beta[c])
      }
    }

  })

  pumpConsts <- list(
    #Model constants
    G = G, B = B, #n = n,                                          #generic
    a  = sum(nx-1), matrix.ind.a = matrix.ind.a,                  #combinations
    Nw = ni, w = nb.nimble$adj, L = sum(ni),                      #spatial
    z_pad = z, C = C,                       #continuous
    #Priori functions constants
    v_mu = 5,                                                     #generic  constants
    v_alfa_kj = 5,                                                #combinations
    a_theta = 0.1, b_theta = 0.1,                                 #spatial constants
    a_phi = 0.1, b_phi = 0.1,                                     #spatial
    a_epsilon = mean(ni)*(0.7^2)*b_epsilon, b_epsilon = b_epsilon, #error
    a_beta = 0.1, b_beta = 0.1                                    #continuous
  )
  pumpData <- list(ci_b = ci_b[1:pumpConsts$G, 1:pumpConsts$B])
  pumpInits <- list(
    #Initial values for the parameters
    mu = 1,                                                       #generic constants
    alfa = rep(0,pumpConsts$a),                                   #combinations
    theta = rep(1, pumpConsts$G),                                 #spatial constants
    phi = matrix(0, pumpConsts$G, pumpConsts$a),                  #spatial
    epsilon = matrix(0, pumpConsts$G, pumpConsts$B),              #error
    beta = matrix(1, pumpConsts$G, pumpConsts$C)                  #continuous
  )

  return(list(Code = pumpCode,
              Consts = pumpConsts,
              Data = pumpData,
              Inits = pumpInits))
}

#Caso contínuo multiple spatial beta####
mult_nimble_beta_spatial <- function(G, B, nx, matrix.ind.a, nb.nimble, ni, z, C, b_epsilon){
  pumpCode <- nimbleCode({
    #Calculation of alfa_b for each combination####
    for (j in 1:B){
      for(k in 1:a){
        sup_alfa[j,k] <- alfa[k] * matrix.ind.a[j,k]
      }

      alfa_b[j] <- sum(sup_alfa[j,1:a])
    }

    #Calculation of phi_b for each combination#####
    for(i in 1:G){
      for (j in 1:B){
        for(k in 1:a){
          sup_phi[i,j,k] <- phi[i,k] * matrix.ind.a[j,k]
        }

        phi_b[i,j] <- sum(sup_phi[i,j,1:a])
      }
    }

    #Sum of the ponderated betas####
    for(i in 1:G){
      for (j in 1:B){
        for (c in 1:C){
          sup_beta[i,j,c] <- beta[i,c]*z_pad[(i+G*(c-1)),j]
        }
        beta_soma[i,j] <- sum(sup_beta[i,j,1:C])
      }
    }

    #Modeling lambda and ci_b#####
    for(i in 1:G){
      for(j in 1:B){
        log(lambda[i,j]) <-  mu +
          alfa_b[j] +
          theta[i] +
          phi_b[i,j] +
          beta_soma[i,j] +
          epsilon[i,j]# + n
        ci_b[i,j] ~ dpois(lambda[i,j])
      }
    }

    #Priori functions####

    for(k in 1:a){
      phi[1:G,k] ~ dcar_normal(adj = w[1:L],
                               num = Nw[1:G],
                               tau = tau_phi[k],
                               zero_mean = 1)
      tau_phi[k] ~ dgamma(a_phi,b_phi)
    }


    theta[1:G] ~ dcar_normal(adj = w[1:L],
                             num = Nw[1:G],
                             tau = tau_theta,
                             zero_mean = 1)
    tau_theta ~ dgamma(a_theta,b_theta)


    for(k in 1:a){
      alfa[k] ~ dnorm(0, var = v_alfa_kj)
    }


    mu ~ dnorm(0, var = v_mu)


    for(i in 1:G){
      for(j in 1:B){
        epsilon[i,j] ~ dnorm(0, tau = tau_e)
      }
    }
    tau_e ~ dgamma(a_epsilon, b_epsilon)

    for(c in 1:C){
      beta[1:G, c] ~ dcar_normal(adj = w[1:L],
                                 num = Nw[1:G],
                                 tau = tau_beta[c],
                                 zero_mean = 1)
      tau_beta[c] ~ dgamma(a_beta, b_beta)
    }

  })

  pumpConsts <- list(
    #Model constants
    G = G, B = B,# n = n,                                         #generic
    a  = sum(nx-1), matrix.ind.a = matrix.ind.a,                  #combinations
    Nw = ni, w = nb.nimble$adj, L = sum(ni),                      #spatial
    z_pad = z, C = C,                                            #continuous
    #Priori functions constants
    v_mu = 5,                                                     #generic  constants
    v_alfa_kj = 5,                                                #combinations
    a_theta = 0.1, b_theta = 0.1,                                 #spatial constants
    a_phi = 0.1, b_phi = 0.1,                                     #spatial
    a_epsilon = mean(ni)*(0.7^2)*b_epsilon, b_epsilon = b_epsilon, #error
    a_beta = 0.1, b_beta = 0.1                                    #continuous
  )
  pumpData <- list(ci_b = ci_b[1:pumpConsts$G, 1:pumpConsts$B])
  pumpInits <- list(
    #Initial values for the parameters
    mu = 1,                                                       #generic constants
    alfa = rep(0,pumpConsts$a),                                   #combinations
    theta = rep(1, pumpConsts$G),                                 #spatial constants
    phi = matrix(0, pumpConsts$G, pumpConsts$a),                  #spatial
    epsilon = matrix(0, pumpConsts$G, pumpConsts$B),              #error
    beta = matrix(1, pumpConsts$G, pumpConsts$C)                  #continuous
  )

  return(list(Code = pumpCode,
              Consts = pumpConsts,
              Data = pumpData,
              Inits = pumpInits))
}
