mkdir <- function(path = NULL){
  if(!is.null(path)&(path!="")){
    for (ndir in 1:length(path)){
      dirtemp=path[ndir]
      dirtemp=gsub("//","/",dirtemp)
      a=as.matrix(strsplit(dirtemp,split = "/")[[1]])
      b=paste(a[1],"/",sep="")
      if(!dir.exists(b)){dir.create(b)}
      if (length(a)>1){
        for (i in 2:length(a)){
          if(nchar(a[i])>0){
            b=paste(b,a[i],"/",sep="")
          }
          if(!dir.exists(b)){dir.create(b)}
        }
      }
    }
  }
}
