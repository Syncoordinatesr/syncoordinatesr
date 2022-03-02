###Função para mu
muf <- function (x,sumeta,vmu,ci_b){
  x*ci_b-((1/(2*vmu))*x^2)-exp(x)*sumeta
}

mufprima <- function(x,sumeta,vmu,ci_b){
  ci_b-(x/vmu)-exp(x)*sumeta
}
