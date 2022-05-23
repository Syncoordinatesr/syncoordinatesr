#' @title  A function that generates useful objects for the package from the database
#'
#' @name  prepare_data
#'
#' @description  An auxiliary function to generate useful objects for the order functions in this package
#' Function \code{prepare_data} has the goal of receiving the database of the user to generate important variables that will be used in the MCMC.
#' and in the end to generate the synthetic coordinates.
#' In the input, the function receives the parameters: \code{dataset}, \code{coord}, \code{grid}.
#'
#' @param   dataset   A data frame with all the information except the coordinates
#' @param   coord   An object with two columns indicating the latitude and longitude respectively of the elements in the dataset
#' @param   grid  The grid represents the quantities of divisions that will be made in the location. Bigger the grid, closer the synthetic coordinates are to the real coordinates. With a default result of (grid = 10)
#'
#' @return  A list containing useful objects to the mcmc function.
#'
#' @references
#' NUNES, Letícia. Métodos de Simulação de Dados Geográficos Sintéticos Para Bases Confidenciais. *Dissertação de Mestrado*, [s. l.], 2018.
#' Disponível em: \url:{http//est.ufmg.br/portal/arquivos/mestrado/dissertacoes/dissertacao_Leticia_Silva_Nunes.pdf}. Acesso em: 2 mar. 2022.
#'
#' @import spdep

prepare_data <- function(dataset, coord, grid = 10){

  # Fazer um tratamento melhor dos erros
  if (!is.data.frame(dataset))
    stop("The dataset must be an object of 'data.frame' class")
  if (!(nrow(dataset) == nrow(coord)))
    stop("Both objects: 'dataset' and 'coord' must have the same numer of lines (elements)")
  #if (!(integer(grid))) stop("The grid should be an integer value")

  n = dim(dataset)[1] ; p = dim(dataset)[2]
  vx = list(1); vx = c(vx,2:p)
  for(i in 1:p){
    vx[[i]] = sort(unique(dataset[,i]))
  }

  nx = numeric(p)
  for(i in 1:p){
    nx[i] = length(vx[[i]])
  }

  B = prod(nx)
  b = 1:B

  b.matrix = matrix(0,B,p+1)
  b.matrix[,1] = b
  b.matrix[,2] = rep(vx[[1]], times=1, each=prod(nx[2:p]))
  for(i in 2:p){
    b.matrix[,i+1] = rep(vx[[i]], times=prod(nx[1:(i-1)]),
                         each=prod(nx[(i+1):p]))
  }

  G = grid*grid

  # Pegando os limites inferiores e superiores da latitude e longitude
  xmax = max(coord[,1])
  xmin = min(coord[,1])
  ymax = max(coord[,2])
  ymin = min(coord[,2])

  # Dividindo em espaços homogêneos de acordo com o 'grid'
  dxvec = (xmax-xmin)/grid
  dyvec = (ymax-ymin)/grid

  # Criando os vetores de latitude e longitude partidas, e criando as células
  lonvec = seq(xmin,xmax,dxvec)
  latvec = seq(ymin,ymax,dyvec)
  xlon = findInterval(coord[,1], lonvec, all.inside = T)
  ylat = findInterval(coord[,2], latvec, all.inside = T)
  celula = matrix(cbind(xlon,ylat), nrow=n, ncol=2)

  cel = numeric(n)
  for (i in 1:n){
    cel[i] = celula[i,1] + (celula[i,2]-1)*grid
  }

  dados.ord = dataset[do.call(order, dataset), ]
  u = unique(dados.ord)

  comb = numeric(n)
  for(i in 1:n){
    for (j in 1:B){
      if (sum(dados[i,1:p]== u[j,])==p)
        comb[i]=j
    }
  }

  c = cbind(cel,comb)

  ci_b = matrix(0,G,B)

  for(i in 1:dim(c)[1]){
    ci_b[c[i,1],c[i,2]] = ci_b[c[i,1],c[i,2]] + 1
  }

  # Criando o fator de influência em células vizinhas
  neigh = cell2nb((grid),(grid), type="queen")
  ni = card(neigh)
  W = nb2mat(neigh, style="B")

  # list of the elements in alpha used in each combination
  ind.a = list(1); ind.a = c(ind.a,2:B)
  for(i in b){
    ind.a[[i]] = numeric(0)
    j=1
    if(any(b.matrix[i,j+1]==vx[[j]][-1])){
      ind.a[[i]] = c(ind.a[[i]], which(b.matrix[i,j+1]==vx[[j]][-1]))
    }
    for(j in 2:p){
      if(any(b.matrix[i,j+1]==vx[[j]][-1])){
        ind.a[[i]] = c(ind.a[[i]], sum(nx[1:(j-1)]-1) +
                         which(b.matrix[i,j+1]==vx[[j]][-1]))
      }
    }
  }

  # list with the combinations that have each alpha
  sub.a = list(1); sub.a = c(sub.a,2:sum(nx-1))
  for(i in 1:sum(nx-1)){
    sub.a[[i]] = numeric(0)
    for(j in b){
      if(i %in% ind.a[[j]])
        sub.a[[i]] = c(sub.a[[i]],j)
    }
  }

  Z <- matrix(,n,1)
  sigma.z=10

  #Falta Z barra

  for (i in 1:n){
    if(comb[i]==1){
      Z[i] <- rnorm(1,zbarra1,sigma.z)}

    if(comb[i]==2){
      Z[i] <- rnorm(1,zbarra2,sigma.z)}

    if(comb[i]==3){
      Z[i] <- rnorm(1,zbarra3,sigma.z)}

    if(comb[i]==4){
      Z[i] <- rnorm(1,zbarra4,sigma.z)}

    if(comb[i]==5){
      Z[i] <- rnorm(1,zbarra5,sigma.z)}

    if(comb[i]==6){
      Z[i] <- rnorm(1,zbarra6,sigma.z)}

    if(comb[i]==7){
      Z[i] <- rnorm(1,zbarra7,sigma.z)}

    if(comb[i]==8){
      Z[i] <- rnorm(1,zbarra8,sigma.z)}

    if(comb[i]==9){
      Z[i] <- rnorm(1,zbarra9,sigma.z)}

    if(comb[i]==10){
      Z[i] <- rnorm(1,zbarra10,sigma.z)}

    if(comb[i]==11){
      Z[i] <- rnorm(1,zbarra11,sigma.z)}

    if(comb[i]==12){
      Z[i] <- rnorm(1,zbarra12,sigma.z)}

  }#Valor deveria depender do valor de b

  for(i in 1:n){
    z.pad[i] <- (Z[i]/mean(Z))/sd(Z)
  }

  return(list(n=n, p=p, vx=vx, nx=nx, B=B, b=b, G=G,
              latvec=latvec, lonvec=lonvec, comb=comb, ci_b=ci_b, ni=ni,
              ind.a=ind.a, sub.a=sub.a, W=W, Z=Z, z.pad=z.pad))
}
