tabulate.char<-function(char,sort=F,decreasing = F){
  #char: a vector of characters
  #sort: logical, sort the output or not

  char=as.matrix(as.character(char))
  unique_char=unique(char)
  if(sum(is.na(unique_char))){
    unique_char=unique_char[!is.na(unique_char)]
    temp=c(N.A=sum(is.na(char)))
    char=char[!is.na(char)]
  }else{
    temp=NULL
  }

  if(length(unique(char))>1){
    char_tabulate=matrix(1,length(unique_char),1,dimnames = list(unique_char,"number"))
    for(i in unique_char){char_tabulate[i,]=sum(char==i)};
    if(sort){
      char_tabulate=sortby(char_tabulate,sort_vector = char_tabulate,decreasing = decreasing)
    }
    if(!is.null(temp)){
      char_tabulate=rbind(char_tabulate,N.A=temp[1])
    }

    colnames(char_tabulate)="number"

  }else{
    if(!is.null(temp)){
      char_tabulate=c(length(char),temp[1]);names(char_tabulate)=c(unique_char,"N.A")
    }else{
      char_tabulate=length(char);names(char_tabulate)=unique_char
    }
  }

  return(char_tabulate)


}
