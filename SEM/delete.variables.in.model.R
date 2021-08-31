delete.variables.in.model<-function(fit=NULL,delete.pars=NULL,delete.paths=NULL,export.model.table=FALSE){
  #delete.paths=list(lhs=c("process"),rhs=c("SF","LF"),op="~")
  par.table <- partable(fit)
  par.table=par.table[par.table$op!=":=",]
  par.table$expression.nolable=with(data = par.table,paste0(lhs,op,rhs))

  no.lable=par.table$label==""
  if(sum(no.lable)){par.table$label[no.lable]=paste("par",1:sum(no.lable),sep="")}

  judge1=par.table$rhs%in%delete.pars | par.table$lhs%in%delete.pars

  if(!is.null(delete.paths)){
    delete.paths=as.data.frame(cbind(lhs=delete.paths$lhs,op=delete.paths$op,rhs=delete.paths$rhs))
    for(i in 1:nrow(delete.paths)){
      if(length(unique(delete.paths[i,c("lhs","rhs")]))==2){
        judge=judge1|(par.table$lhs==delete.paths[i,"lhs"]&par.table$rhs==delete.paths[i,"rhs"]&par.table$op==delete.paths[i,"op"])
      }
    }
  }else{
    judge=judge1
  }
  if(sum(judge)){
    print(paste0(sum(judge)," path",if(sum(judge)>1){"s were "}else{" was "},"deleted. < ",
                 paste(par.table$expression.nolable[judge],collapse = " & ")," >" ))
    par.table=par.table[!judge,]
  }else{
    warning("no paths deleted")
  }



  model=paste0(sort(paste(par.table$lhs,par.table$op,par.table$rhs,sep="")),collapse ="; ")
  if(export.model.table){
    model.table=par.table[,c("lhs","op","rhs","label")]
    model.table$experession=paste(par.table$lhs,par.table$op,par.table$label,"*",par.table$rhs,sep="")
    model.table$experession.nolable=paste(par.table$lhs,par.table$op,par.table$rhs,sep="")
    rownames(model.table)=model.table$experession.nolable
    return(list(model=model,model.table=model.table))
  }else{
    return(model=model)
  }
}

