char.align<-function(x,align=c("none","left","right"),digits=3){
  #(l=char.align(c(rep(NA,10),runif(100)/10),align = "left"))
  #(r=char.align(c(rep(NA,10),runif(100)/10),align = "right"))
  #paste(l,r,sep=" - ")
  align.X=vector(length = length(x));
  ISNA.X=is.na(x)
  if(sum(ISNA.X)){align.X[ISNA.X]="NA"}
  
  if(sum(!ISNA.X)){x=as.vector(x[!ISNA.X])}
  if(is.numeric(x)){x=round(x,digits)}
  if(!is.character(x)){x=as.character(x)}
  
  align.X[!ISNA.X]=x
  
  if(align[1]=="right"){
    temp=max(nchar(align.X))-nchar(align.X)
    align.char=apply(as.matrix(temp),1,FUN=function(x){return(if(x==0){""}else{paste0(rep(" ",x),collapse = "")})})
    align.X=paste0(align.char,align.X)
  }else if(align[1]=="left"){
    temp=max(nchar(align.X))-nchar(align.X)
    align.char=apply(as.matrix(temp),1,FUN=function(x){return(if(x==0){""}else{paste0(rep(" ",x),collapse = "")})})
    align.X=paste0(align.X,align.char)
  }
  return(align.X)
}

           