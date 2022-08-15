#' @title  A function that generates useful objects for the package from the database
#'
#' @name  prepare_data
#'
#' @description  An auxiliary function to generate useful objects for the other functions in this package.
#' Function \code{prepare_data} has the goal of receiving the database of the user to generate important variables that will be used in the MCMC.
#' And in the end to generate the synthetic coordinates.
#' In the input, the function receives the parameters: \code{dataset}, \code{coord}, \code{grid},  \code{continuous}.
#'
#' @param   dataset   A data frame with all the information except the coordinates
#' @param   coord   An object with two columns indicating the latitude and longitude respectively of the elements in the dataset
#' @param   grid  The grid represents the quantities of divisions that will be made in the location. Bigger the grid, closer the synthetic coordinates are to the real coordinates. With a default result of (grid = 10)
#' @param   continuous  An object that indicates which columns in the dataset correspond to continuous variables. The default is FALSE which means that there is none continuous variable. (Still not adapted for cases with more than one continuous variable)
#'
#' @return  A list containing useful objects to the mcmc function.
#'
#' @references
#' NUNES, Letícia. Métodos de Simulação de Dados Geográficos Sintéticos Para Bases Confidenciais. *Dissertação de Mestrado*, [s. l.], 2018.
#' Disponível em: \url:{http//est.ufmg.br/portal/arquivos/mestrado/dissertacoes/dissertacao_Leticia_Silva_Nunes.pdf}. Acesso em: 2 mar. 2022.
#'
#' @import spdep

prepare_data <- function(dataset, coord, grid = 10, continuous = FALSE){

  if (!is.data.frame(dataset))
    stop("The dataset must be an object of 'data.frame' class")
  if (!(nrow(dataset) == nrow(coord)))
    stop("Both objects: 'dataset' and 'coord' must have the same numer of lines (elements)")
  #if (!(integer(grid))) stop("The grid should be an integer value")
  if (continuous == TRUE)
    stop("The continuous object must indicates which columns in the dataset correspond to continuous variables")


  n = dim(dataset)[1]
  n_coord = dim(coord)[1]
  if(n != n_coord)
    stop("Both objects 'dataset' and 'coord' must have the same amount of elements")


  if(continuous != FALSE){

    vZ = length(continuous)
    Z = matrix(NA, n, vZ)

    for(i in 1:vZ){
      col_Z = continuous[i]
      Z[, i] = dataset[, col_Z]
    }
    for(i in 1:vZ){
      col_Z = continuous[i]
      dataset = dataset[, -col_Z]
    }
  }

  p = dim(dataset)[2]

  vx = list(1); vx = c(vx,2:p)
  for(i in 1:p){
    vx[[i]] = sort(unique(dataset[,i]))
  }

  nx = numeric(p)
  for(i in 1:p){
    nx[i] = length(vx[[i]])
  }

  #B = prod(nx)
  #b = 1:B

  #b.matrix = matrix(0,B,p+1)
  #b.matrix[,1] = b
  #b.matrix[,2] = rep(vx[[1]], times=1, each=prod(nx[2:p]))
  #for(i in 2:p){
  #  b.matrix[,i+1] = rep(vx[[i]], times=prod(nx[1:(i-1)]),
  #                       each=prod(nx[(i+1):p]))
  #}

  G = grid*grid

  # Lower and upper limits of latitude and longitude
  xmax = max(coord[,1])
  xmin = min(coord[,1])
  ymax = max(coord[,2])
  ymin = min(coord[,2])

  # Dividing in equal spaces according to the 'grid' value
  dxvec = (xmax-xmin)/grid
  dyvec = (ymax-ymin)/grid

  # Vectors of latitude and longitude divided
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
  b.matrix = unique(dados.ord)

  B = nrow(b.matrix)
  b = 1:B
  b.matrix = cbind(b, b.matrix)

  comb = numeric(n)
  for(i in 1:n){
    for (j in 1:B){
      if (sum(dataset[i,1:p]== b.matrix[j,-1])==p)
        comb[i]=j
    }
  }

  c = cbind(cel,comb)

  ci_b = matrix(0,G,B)

  for(i in 1:dim(c)[1]){
    ci_b[c[i,1],c[i,2]] = ci_b[c[i,1],c[i,2]] + 1
  }

  # The influence factor in neighboring cells
  neigh = cell2nb((grid),(grid), type="queen")
  ni = card(neigh)
  W = nb2mat(neigh, style="B")

  # List of the elements in alpha used in each combination
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


  if(continuous != FALSE){

    z.bar = z.pad = array(NA, c(G, B, vZ))

    for(d in 1:vZ){
      for(i in 1:G){
        for(j in 1:B){
            z.bar[i, j, d] = ifelse(sum(cel == i & comb == j) == 0,
                                    NA,
                                    mean(Z[cel == i & comb == j, d], na.rm = TRUE))
        }
      }
    }

    for(d in 1:vZ){
      for(b in 1:B){
        z.pad[, b, d] = (z.bar[, b, d]-mean(z.bar[, b, d], na.rm = TRUE))/sd(z.bar[, b, d], na.rm = TRUE)
      }
    }
  }

  if(continuous != FALSE){
    return(list(n=n, p=p, vx=vx, nx=nx, B=B, b=b, G=G,
                latvec=latvec, lonvec=lonvec, comb=comb, ci_b=ci_b, ni=ni,
                ind.a=ind.a, sub.a=sub.a, W=W, Z=Z, vZ=vZ, z.pad=z.pad))
  } else{
    return(list(n=n, p=p, vx=vx, nx=nx, B=B, b=b, G=G,
                latvec=latvec, lonvec=lonvec, comb=comb, ci_b=ci_b, ni=ni,
                ind.a=ind.a, sub.a=sub.a, W=W))
  }
}
