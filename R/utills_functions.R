# utills_functions

## Função para beta_i_ICAR
logit <- function(x){
  y <- exp(x)/(1+exp(x))
  y
}

inv.logit <-function(x){
  y <- log(x/(1-x))
  y
}

betaf_ICAR = function(x, c.f, zib.vec, eta,ni,tau.beta,bar.f){
  bb <- inv.logit(x)
  bb*c.f - sum( exp(bb*zib.vec + eta) ) - (1/2)*ni*tau.beta*(bb-bar.f)^2 - log(x) - log(1-x)
}
##

## Função para beta_i_Normal
logit <- function(x){
  y <- exp(x)/(1+exp(x))
  y
}

inv.logit <-function(x){
  y <- log(x/(1-x))
  y
}

betaf = function(x, c.f, zib.vec, eta, vbeta){
  bb <- inv.logit(x)
  bb*c.f - sum( exp(bb*zib.vec + eta) ) - (1/(2*vbeta))*bb^2 - log(x) - log(1-x)
}
##

## Função para criar as medidas de distâncias
## Usada para avaliar o risco dos modelos
cria.d = function(mat){
  nb = dim(mat)[1]
  d = numeric(nb)
  for(i in 1:nb){
    pos = which(mat==min(mat),arr.ind = T)
    d[i] = mat[pos]
    if(i==nb)
      break
    if(i<(nb-1)){
      mat = mat[-pos[1],]
      mat = mat[,-pos[2]]
    }
    else{
      mat = mat[-pos[1],]
      mat = mat[-pos[2]]
    }
  }
  return(d)
}
##

## Função para epsilon_i
ef = function(x,ci_b,sumeta,tau.e){
  x*ci_b - exp(x)*sumeta - (1/2)*tau.e*x^2
}

efprima = function(x,ci_b,sumeta,tau.e){
  ci_b - exp(x)*sumeta - tau.e*x
}
##

## Function to verify if the points are inside the interval passed by the user (limits)
is_in_area <- function(x, y, xmin, xmax, ymin, ymax){
  is_x_valid = (x >= xmin & x <= xmax)
  is_y_valid = (y >= ymin & y <= ymax)
  return((is_x_valid & is_y_valid))
}
##

## Função para mu
muf <- function (x,sumeta,vmu,ci_b){
  x*ci_b-((1/(2*vmu))*x^2)-exp(x)*sumeta
}

mufprima <- function(x,sumeta,vmu,ci_b){
  ci_b-(x/vmu)-exp(x)*sumeta
}
##

## Função para theta_i
thetaf = function(x, ci_b, sumeta, ni, tau.f, bar.f){
  # x*ci_b - exp(x)*sumeta - (1/2)*ni*tau.f*(x-bar.f)^2
  x*ci_b - exp(x)*sumeta - (1/2)*ni*tau.f*(x-as.vector(bar.f))^2
}

# updated by Thais (Feb/2022) after warning messages about dimension of bar.f

thetafprima = function(x, ci_b, sumeta, ni, tau.f, bar.f){
  # ci_b - exp(x)*sumeta - ni*tau.f*(x-bar.f)
  ci_b - exp(x)*sumeta - ni*tau.f*(x-as.vector(bar.f))
}
##

## Functions to convert dimensions between matrix and vectors
## Thais Paiva, Sep 2011

# vec2mat
# Function to recover the matrix index
vec2mat = function(i,size){
  row = ceiling(i/(size))
  col = i - (row-1)*(size)
  matrix(c(col,row),length(i),2,byrow=F)
}

# mat2vec
# Function to recover the vector index
mat2vec = function(r,c,size){
  matrix((r-1)*(size) + c,length(r),1)
}
##

