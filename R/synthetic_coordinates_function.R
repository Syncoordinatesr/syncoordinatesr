#' @title  Generation of synthetic coordinates
#'
#' @name  syncoordinates
#'
#' @description  A function that generates synthetic coordinates.
#' Function \code{syncoordinates} receives the database, the parameter lambda and the number of synthetic data the user desires.
#' And the function returns the synthetic databases containing the synthetic coordinates.
#'
#' @param  dataset  A data frame with all the information except the coordinates
#' @param  coord  An object with two columns indicating the longitude and latitude respectively of the elements in the dataset
#' @param  grid  The grid represents the quantities of divisions that will be made in the location. Bigger the grid, closer the synthetic coordinates are to the real coordinates. With a default result of (grid = 10)
#' @param  continuous  An object that indicates which columns in the dataset correspond to continuous variables. The default is FALSE which means that there is none continuous variable. (Still not adapted for cases with more than one continuous variable)
#' @param  list_mcmc  Output of the mcmc function
#' @param  n.syn  Number of synthetic database that will be returned
#'
#' @return The return will depend on the argument continuous, if '\code{continuous} = FALSE' the function will return an object of \code{data.frame} class containing all new synthetic coordinates, but if '\code{continuous} != FALSE' in addition to the \code{data.frame} with the synthetic coordinates the function will return new synthetic data to each continuous variables indicated in argument \code{continuous}
#'
#' @references
#' NUNES, Letícia. Métodos de Simulação de Dados Geográficos Sintéticos Para Bases Confidenciais. *Dissertação de Mestrado*, [s. l.], 2018.
#' Disponível em: \url:{http//est.ufmg.br/portal/arquivos/mestrado/dissertacoes/dissertacao_Leticia_Silva_Nunes.pdf}. Acesso em: 2 mar. 2022.
#'
#' @examples
#'   syncoordinates(dataset = my_database, lambda = lambda_from_mcmc, n.syn = 3)
#'
#' @export

syncoordinates <- function(dataset, coord, grid = 10, continuous = FALSE, list_mcmc, n.syn = 5){

  saida = prepare_data(dataset, coord, limits = c(), grid, continuous)

  mapply(assign, names(saida), saida, MoreArgs=list(envir = globalenv()))
  mapply(assign, names(list_mcmc), list_mcmc, MoreArgs=list(envir = globalenv()))

  int.syn = floor((S-burn-1)/(n.syn-1))
  list.syn = seq(burn+1, burn+(n.syn-1)*int.syn+1, by=int.syn)
  #syn.data = list(1)
  #syn.data = c(syn.data,2:n.syn)

  m = length(list.syn)

  prob = array(0, dim=c(G,m,B))

  for(k in 1:m){
    for(j in 1:B){
      for(i in 1:G){
        prob[i,k,j] = lambda[i,list.syn[k],j]/sum(lambda[, list.syn[k], j])
      }
    }
  }

  loc = matrix(0, n, m)

  for(k in 1:m){
    for (j in 1:B){
      aux = comb==j
      loc[aux,k] = sample(c(1:G), sum(aux), prob=prob[,k,j], rep=T)
    }
  }

  coord.syn = array(0, c(n,2,m))
  aux1 = array(0, c(n,2,m))

  for(k in 1:m){
    aux1[,,k] = vec2mat(loc[,k], grid)
    for (i in 1:n){
      coord.syn[i,1,k] = runif(1, lonvec[aux1[i,1,k]], lonvec[aux1[i,1,k]+1])
      coord.syn[i,2,k] = runif(1, latvec[aux1[i,2,k]], latvec[aux1[i,2,k]+1])
    }
  }

  if(continuous != FALSE){

    sigma.z = matrix(0, B, vZ)
    mean.z = matrix(0, B, vZ)

    for(i in 1:B){
      for(j in 1:vZ){
        sigma.z[i,j] = sd(z.bar[, i, j], na.rm = T)
        mean.z[i,j] = mean(z.bar[, i, j], na.rm = T)
        #sigma.z[i,j] = sd(z.bar[z.bar[,i,j]>0, i,])
        #mean.z[i,j] = mean(z.bar[z.bar[,i,j]>0, i,])
      }
    }

    z.syn = array(0, c(n, vZ, m))

    for(k in 1:m){
      for(j in 1:vZ){
        for(i in 1:n){
          if(is.na(z.pad[loc[i,k],comb[i],j]) == FALSE){
            z.syn[i,j,k] = rnorm(1,z.pad[loc[i,k],comb[i],j],1) * sigma.z[comb[i],j] + mean.z[comb[i],j]
          }
          #z.syn[i,j,k] = rnorm(1,0,1) * sigma.z[comb[i],j] + mean.z[comb[i],j]
          else{
            z.pad[loc[i,k],comb[i],j] = 0
            z.syn[i,j,k] = rnorm(1,z.pad[loc[i,k],comb[i],j],1) * sigma.z[comb[i],j] + mean.z[comb[i],j]
          }
        }
      }
    }

    #dataset.syn <- array(NA, c(n, p, m))
    #aux.dataset <- dataset
    #
    #for(k in 1:m){
    #   for(j in 1:vZ){
    #      aux.dataset[, continuous[j]] <- z.syn[,j,k]
    #   }
    #   dataset.syn[,,k] <- aux.dataset
    #}

  }

  if(is.numeric(continuous)){
    return(list(coord.syn, z.syn))
  } else{
    return(coord.syn)
  }
}
