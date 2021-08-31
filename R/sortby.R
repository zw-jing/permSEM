sortby<-function(x,sort_vector=1, MARGIN = 2, decreasing = FALSE, na.last = T){
  #sort a matrix (vector) by a single row or column
  #x: input matrix or vector
  #sort_vector, by which cols (colname) or rows (rowname) to sort the matrix

  if (is.vector(x)){x=as.matrix(x)}

  if(length(sort_vector)==1){sort_vector=if(MARGIN == 2){x[,sort_vector]}else{x[sort_vector,]}}
  if (MARGIN==1){sorted=x[,order(sort_vector,decreasing=decreasing)]
  }else{
    sorted=x[order(sort_vector,decreasing=decreasing),]
  }
  if(is.vector(sorted)){sorted=as.matrix(sorted)}
  if(!is.null(rownames(x))&is.null(rownames(sorted))){rownames(sorted)=rownames(x)}
  if(!is.null(colnames(x))&is.null(colnames(sorted))){colnames(sorted)=colnames(x)}
  return(sorted)
}


