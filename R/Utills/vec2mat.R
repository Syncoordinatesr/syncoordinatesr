#### Functions to convert dimensions between matrix and vectors
#### Thais Paiva, Sep 2011

## vec2mat
## function to recover the matrix index

vec2mat = function(i,size){
  row = ceiling(i/(size))
  col = i - (row-1)*(size)
  matrix(c(col,row),length(i),2,byrow=F)
}

## mat2vec
## function to recover the vector index


mat2vec = function(r,c,size){
  matrix((r-1)*(size) + c,length(r),1)
}
