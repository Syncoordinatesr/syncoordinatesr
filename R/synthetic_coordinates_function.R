#' @title  Generation of synthetic coordinates
#'
#' @name  syncoordinates
#'
#' @description  A function that generates synthetic coordinates.
#' Function \code{syncoordinates} receives the database, the parameter lambda and the number of synthetic data the user desires.
#' And the function returns the synthetic databases containing the synthetic coordinates.
#'
#' @param  dataset   A data frame with all the information except the coordinates
#' @param  coord    An object with two columns indicating the latitude and longitude respectively of the elements in the dataset
#' @param  grid    The grid represents the quantities of divisions that will be made in the location. Bigger the grid, closer the synthetic coordinates are to the real coordinates. With a default result of (grid = 10)
#' @param  list_mcmc    Output of the mcmc function
#' @param  n.syn   number of synthetic database that will be returned
#'
#' @return The function returns an object of \code{data.frame} class, containing all new synthetic coordinates
#'
#' @references
#' NUNES, Letícia. Métodos de Simulação de Dados Geográficos Sintéticos Para Bases Confidenciais. *Dissertação de Mestrado*, [s. l.], 2018.
#' Disponível em: \url:{http//est.ufmg.br/portal/arquivos/mestrado/dissertacoes/dissertacao_Leticia_Silva_Nunes.pdf}. Acesso em: 2 mar. 2022.
#'
#' @examples
#'   syncoordinates(dataset = my_database, lambda = lambda_from_mcmc, n.syn = 3)
#'
#' @export

syncoordinates <- function(dataset, coord, grid = 10, list_mcmc, n.syn = 5){

  saida = prepare_data(dataset, coord, grid)

  int.syn = floor((list_mcmc$S-list_mcmc$burn-1)/(n.syn-1))
  list.syn = seq(list_mcmc$burn+1,list_mcmc$burn+(n.syn-1)*int.syn+1,by=int.syn)
  #syn.data = list(1)
  #syn.data = c(syn.data,2:n.syn)

  m=length(list.syn)

  prob=array(0,dim=c(saida$G,m,saida$B))

  for(k in 1:m){
    for(j in 1:saida$B){
      for(i in 1:saida$G){
        prob[i,k,j]=list_mcmc$lambda[i,list.syn[k],j]/sum(list_mcmc$lambda[,list.syn[k],j])
      }
    }
  }

  loc=matrix(0,saida$n,m)

  for(k in 1:m){
    for (j in 1:saida$B){
      aux = saida$comb==j
      loc[aux,k] = sample(c(1:saida$G),sum(aux),prob=prob[,k,j],rep=T)
    }
  }

  sin=array(0,c(saida$n,2,m))
  aux1=array(0,c(saida$n,2,m))

  for(k in 1:m){
    aux1[,,k] = vec2mat(loc[,k],grid)
    for (i in 1:saida$n){
      sin[i,1,k] = runif(1, saida$lonvec[aux1[i,1,k]], saida$lonvec[aux1[i,1,k]+1])
      sin[i,2,k] = runif(1, saida$latvec[aux1[i,2,k]], saida$latvec[aux1[i,2,k]+1] )
    }
  }
  return(sin)
}
