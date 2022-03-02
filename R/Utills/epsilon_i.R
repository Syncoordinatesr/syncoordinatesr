###Função para epsilon_i
ef = function(x,ci_b,sumeta,tau.e){
  x*ci_b - exp(x)*sumeta - (1/2)*tau.e*x^2
}

efprima = function(x,ci_b,sumeta,tau.e){
  ci_b - exp(x)*sumeta - tau.e*x
}
