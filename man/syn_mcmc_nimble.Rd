% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_nimble.R
\name{syn_mcmc_nimble}
\alias{syn_mcmc_nimble}
\title{Application of the Markov Chain Monte Carlo using the nimble package}
\usage{
syn_mcmc_nimble(
  dataset,
  coord,
  limits = c(),
  grid = 10,
  S = 5000,
  burn = 1000,
  continuous = FALSE,
  spatial_beta = FALSE,
  return_parameters = FALSE
)
}
\arguments{
\item{dataset}{A data frame with all the information except the coordinates}

\item{coord}{An object with two columns indicating the latitude and longitude respectively of the elements in the dataset}

\item{limits}{An object that is a vector of the dimensions where will be create the grids passed through the sequence of xmin, xmax, ymin, ymax. The default is create by using the maximum and the minimum of the coords object.}

\item{grid}{The grid represents the quantities of divisions that will be made in the location. Bigger the grid, closer the synthetic coordinates are to the real coordinates. With a default result of (grid = 10)}

\item{S}{Quantities of simulations that will be made. With a default result of (S = 5000)}

\item{burn}{The number of simulations that will be burned to warm-up the \code{mcmc}. With a default result of (\code{burn} = 1000)}

\item{continuous}{Option so the user can warn the function for the presence of continuous variables in the dataset. You must indicate the columns numbers that contain continuous variables in a vector}

\item{spatial_beta}{Option so the user can choose to use a spatial \code{beta} parameter. You must input a vector containing which \code{beta} you want to add the spatiality. If you want to select all beta, you can use the logical argument TRUE}

\item{return_parameters}{Option to return the result of the parameters. With default not returning the parameters}

\item{nchain}{Quantities of chains that will be made. With a default result of (nchain = 1)}
}
\value{
Depending on the \code{return_parameters} parameter, this function can return only the \code{lambda} parameter or all other significant parameters too.
}
\description{
To obtain the synthetic coordinates we are going to use the package nimble to apply the MCMC.
Function "syn_mcmc_nimble" receives our database uses variables that we got in the function \code{link:prepare_data} to do the mcmc with nimble.
By the end of this function we get the parameter \code{lambda} that will be required when creating the synthetic coordinates.
}
\examples{
  syn_mcmc_nimble(dataset = my_database, S = 2500, burn = 500, return_parameters = TRUE)

}
\references{
NUNES, Letícia. Métodos de Simulação de Dados Geográficos Sintéticos Para Bases Confidenciais. *Dissertação de Mestrado*, [s. l.], 2018.
Disponível em: \url:{http//est.ufmg.br/portal/arquivos/mestrado/dissertacoes/dissertacao_Leticia_Silva_Nunes.pdf}. Acesso em: 2 mar. 2022.
}
