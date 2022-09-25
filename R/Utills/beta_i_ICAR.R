# Função para beta_i
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
