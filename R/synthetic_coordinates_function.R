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
#' @param  limits  Object that is a vector of the dimensions where the grid cells will be created, specified by the sequence xmin, xmax, ymin, ymax. The default is to use the maximum and minimum of the coords object.
#' @param  grid  The grid represents the quantities of divisions that will be made in the location. Bigger the grid, closer the synthetic coordinates are to the real coordinates. With a default result of (grid = 10)
#' @param  continuous  An object that indicates which columns in the dataset correspond to continuous variables. The default is FALSE which means that there is none continuous variable. (Still not adapted for cases with more than one continuous variable)
#' @param  restricted_area  Logical argument indicating whether there are restricted areas where synthetic geographic coordinates should not be generated. The default is FALSE, meaning there are no restricted areas.
#' @param  coord_restricted_area  Object with two columns indicating the longitude and latitude of the points that form the restricted areas. By default, for the object to be considered a polygon, the first longitude and latitude points must be identical to the last longitude and latitude points.
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
#'    syncoordinates(dataset = my_database , coord = my_coords, grid = 10, continuous = FALSE, restricted_area = FALSE, list_mcmc = my_mcmc, n.syn = 5)
#' @export

syncoordinates <- function(dataset, coord, limits = c(), grid = 10,
                           continuous = FALSE,
                           restricted_area = FALSE,
                           coord_restricted_area,
                           list_mcmc, n.syn = 5){


  saida = prepare_data(dataset, coord, limits, grid, continuous)
  mapply(assign, names(saida), saida, MoreArgs=list(envir = globalenv()))
  mapply(assign, names(list_mcmc), list_mcmc, MoreArgs=list(envir = globalenv()))


  int.syn = floor((S-burn-1)/(n.syn-1))

  list.syn = seq(burn+1, burn+(n.syn-1)*int.syn+1, by=int.syn)


  cat("Dimensões de lambda:\n")
  print(dim(lambda))
  cat("Lista de list.syn:\n")
  print(list.syn)

  m = length(list.syn)
  prob = array(0, dim=c(G,m,B))

  if(restricted_area){
  coord_restricted_area = as.matrix(coord_restricted_area)

  initial_point = FALSE
  polygons = c()
  for (i in 1:dim(coord_restricted_area)[1]) {
    if (class(initial_point) == "logical") {
      initial_point = coord_restricted_area[i,]
      sub_polygons =  c(c(coord_restricted_area[i,1], coord_restricted_area[i,2]))
    } else if ((coord_restricted_area[i,1] == initial_point[1]) & (coord_restricted_area[i,2] == initial_point[2])) {
      sub_polygons = c(sub_polygons, c(coord_restricted_area[i,1], coord_restricted_area[i,2]))
      sub_polygons = matrix(sub_polygons, ncol=2, byrow=TRUE)
      sub_polygons = st_polygon(list(sub_polygons))
      polygons = c(polygons, sub_polygons)
      initial_point = FALSE
      sub_polygons = c()
    } else {
      sub_polygons = c(sub_polygons, c(coord_restricted_area[i,1], coord_restricted_area[i,2]))
    }
  }

  n_polygons = length(polygons)
  cat("Warning: Number of polygons informed in the restricted area = ",
      n_polygons, "\n")
  restricted_area_sf = st_multipolygon(list(polygons))
  }

  if(restricted_area != FALSE){
    lon_vector = lonvec
    lat_vector = latvec
    a = matrix(0, grid, grid)
    size_lon_vector = length(lon_vector) - 1
    size_lat_vector = length(lat_vector) - 1
    intersect_list = c()
    grid_list = c()

    for(j in 1:size_lat_vector){
      for(i in 1:size_lon_vector){
        coords = rbind(c(lon_vector[i], lat_vector[j]),
                       c(lon_vector[i+1], lat_vector[j]),
                       c(lon_vector[i+1], lat_vector[j+1]),
                       c(lon_vector[i], lat_vector[j+1]),
                       c(lon_vector[i], lat_vector[j]))
        grid_cell_sf = st_polygon(list(coords))
        grid_list = c(grid_list, grid_cell_sf)
      }
    }
    grid_polygons = st_multipolygon(list(grid_list))
    for(j in 1:size_lat_vector){
      for(i in 1:size_lon_vector){
        intersect_sum = 0
        cells = st_polygon(grid_polygons[[1]][(j-1)*grid+i])
        for(ra in restricted_area_sf[[1]]){
          ra = st_polygon(list(ra))
          intersect_aux = st_intersection(cells, ra)
          intersect_list = c(intersect_list, intersect_aux)
          intersect_sum = intersect_sum + st_area(intersect_aux)
          intersect_sum = intersect_sum + (st_area(st_intersection(cells, ra)))
        }
        a[j, i] = 1 - intersect_sum / st_area(cells)
      }
    }
    a_vec = c(t(a))
  }else{
    a_vec = rep(1, (grid*grid))
  }

  for(k in 1:m){
    for(j in 1:B){
      sum_lambda = sum(lambda[, list.syn[k], j])
      for(i in 1:G){
        prob[i,k,j] = ( lambda[i,list.syn[k],j] * a_vec[i] ) / ( sum_lambda * a_vec[i] )
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
      flag = TRUE
      while(flag){
        longitude_candidate = runif(1, lonvec[aux1[i,1,k]], lonvec[aux1[i,1,k]+1])
        latitude_candidate = runif(1, latvec[aux1[i,2,k]], latvec[aux1[i,2,k]+1])

        ponto = st_point(c(longitude_candidate, latitude_candidate))

        if(restricted_area){
          arr = array(0, length(restricted_area_sf[[1]]))

          idx = 1

          for(ra in restricted_area_sf[[1]]){
              ra = st_polygon(list(ra))
              resultado = st_intersects(ponto, ra)
              arr[idx] = !is.na(resultado[1] == 1)
              idx = idx + 1
          }

          flag = sum(arr)
        }else{
          flag = FALSE
        }
      }
      coord.syn[i,1,k] = longitude_candidate
      coord.syn[i,2,k] = latitude_candidate
    }
  }


  if(is.numeric(continuous)){

    sigma.z = matrix(0, B, vZ)
    mean.z = matrix(0, B, vZ)

    for(i in 1:B){
      for(j in 1:vZ){
        sigma.z[i,j] = sd(z.bar[, i, j], na.rm = T)
        mean.z[i,j] = mean(z.bar[, i, j], na.rm = T)
      }
    }

    z.syn = array(0, c(n, vZ, m))

    for(k in 1:m){
      for(j in 1:vZ){
        for(i in 1:n){
          if(is.na(z.pad[loc[i,k],comb[i],j]) == FALSE){
            z.syn[i,j,k] = rnorm(1,z.pad[loc[i,k],comb[i],j],1) * sigma.z[comb[i],j] + mean.z[comb[i],j]
          }
          else{
            z.pad[loc[i,k],comb[i],j] = 0
            z.syn[i,j,k] = rnorm(1,z.pad[loc[i,k],comb[i],j],1) * sigma.z[comb[i],j] + mean.z[comb[i],j]
          }
        }
      }
    }

  }

  if(is.numeric(continuous)){
    return(list(coord.syn = coord.syn, z.syn =  z.syn))
  } else{
    return(list(coord.syn))
  }
}
