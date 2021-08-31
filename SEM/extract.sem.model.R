extract.sem.model<-function(fit=NULL,export.model.table=TRUE){
  par.tab <- parTable(fit)
  no.lable=par.tab$label==""
  if(sum(no.lable)){par.tab$label[no.lable]=paste("par",1:sum(no.lable),sep="")}
  model.with.self_def=paste0(sort(paste(par.tab$lhs,par.tab$op,par.tab$label,"*",par.tab$rhs,sep="")),collapse ="; ")
  par.tab=par.tab[par.tab$op%in%c("~","~~"),]
  model=paste0(sort(paste(par.tab$lhs,par.tab$op,par.tab$label,"*",par.tab$rhs,sep="")),collapse ="; ")
  model.nolable=paste0(sort(paste(par.tab$lhs,par.tab$op,par.tab$rhs,sep="")),collapse ="; ")
  if(export.model.table){

    model.table=par.tab[,c("lhs","op","rhs","label")]
    model.table$experession=paste(par.tab$lhs,par.tab$op,par.tab$label,"*",par.tab$rhs,sep="")
    model.table$experession.nolable=paste(par.tab$lhs,par.tab$op,par.tab$rhs,sep="")
    rownames(model.table)=model.table$experession.nolable
    return(list(model=model,model.with.self_def=model.with.self_def,model.table=model.table,model.nolable=model.nolable))
  }else{
    return(list(model=model,model.with.self_def=model.with.self_def,model.nolable=model.nolable))
  }
}


