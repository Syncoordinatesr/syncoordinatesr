###Função para theta_i
thetaf = function(x, ci_b, sumeta, ni, tau.f, bar.f){
  # x*ci_b - exp(x)*sumeta - (1/2)*ni*tau.f*(x-bar.f)^2
  x*ci_b - exp(x)*sumeta - (1/2)*ni*tau.f*(x-as.vector(bar.f))^2
}

## updated by Thais (Feb/2022) after warning messages about dimension of bar.f

thetafprima = function(x, ci_b, sumeta, ni, tau.f, bar.f){
  # ci_b - exp(x)*sumeta - ni*tau.f*(x-bar.f)
  ci_b - exp(x)*sumeta - ni*tau.f*(x-as.vector(bar.f))
}
