# 'Syncoordinatesr' package

The *syncoordinatesr* is a package useful to generate synthetic coordinates for your database. 
Inspired by the research article of Thaís Paiva: <https://onlinelibrary.wiley.com/doi/10.1002/sim.6078>. And the masters dissertation of Letícia Nunes: <http://est.ufmg.br/portal/arquivos/mestrado/dissertacoes/dissertacao_Leticia_Silva_Nunes.pdf>.
The package contains functions able to generate synthetic coordinates using the method explained in both articles. This method simulates synthetic data, that are generated from specified probability distributions with the intention to be as close as possible of the original data, but without spreading the original information and don't compromising the final analysis. 

## Installing the package

The package, in this moment, can't be found on CRAN and your installation has to
be done via github. To do this you must have already installed the package **devtools**
so you can use the function *install_github* as exhibited below:

```R
devtools::install_github("Syncoordinatesr/syncoordinatesr")
```

## Error in installation

In case of error, please contact me on my e-mail: leomgal20@hotmail.com
