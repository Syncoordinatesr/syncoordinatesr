########## Função para criar as medidas de distâncias
#### Usada para avaliar o risco dos modelos

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
