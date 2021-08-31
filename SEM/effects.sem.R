effects.sem<-function(rh=NULL,lh=NULL,fit,model=NULL){
  #effects.model(rh=c("D","I","PR"),lh=c("Multi.stibility","M_Wellbeing"),edges,model=MD)
  par.tab_model=extract.sem.model(fit)
  par.tab=par.tab_model$model.table
  edges=par.tab[par.tab$op=="~",]
  
  if(is.null(rh)){
    if(length(unique(edges$lhs))==1){
      rh=unique(edges$rhs)
    }else{
      rh=rownames(sortby(tabulate.char(edges$rhs),sort_vector = 1,decreasing = T))
    }

  }
  if(is.null(lh)){
    if(length(unique(edges$lhs))==1){
      lh=unique(edges$lhs)
    }else{
      lh=rownames(sortby(tabulate.char(edges$lhs),1))
    }

  }
    
  j1=rh%in%lh
  if(sum(j1)>0&sum(j1)<length(rh)){rh=c(rh[!j1],rh[j1])}
  j2=lh%in%rh
  if(sum(j2)>0&sum(j2)<length(lh)){lh=c(lh[j2],lh[!j2])}
  

  
  
  if(is.null(model)){
    model=par.tab_model$model
  }
  
  mat=matrix(NA,length(rh),length(lh),dimnames = list(rh,lh))
  total.eff=matrix("",length(rh),length(lh),dimnames = list(rh,lh))
  direct.eff=indirect.eff=total.eff
  total.eff.id=indirect.eff.id=direct.eff.id=total.eff
  for(i in 1:length(rh)){
    #if(rh_i==rh[3]&lh_j==lh[1]){break}
    rh_i=rh[i]

    for(j in 1:length(lh)){
      lh_j=lh[j]

      effects=SEM.total.effect.expression(rh=rh_i,lh=lh_j,edges)
      total.eff.id[rh_i,lh_j]=paste("total_rh",i,"lh",j,sep="")
      total.eff[rh_i,lh_j]=if(!is.null(effects$total.effect.elem)){paste(total.eff.id[rh_i,lh_j],":=",paste0(effects$total.effect.elem,collapse ="+"),sep="")}else{""}
      
      indirect.eff.id[rh_i,lh_j]=paste("indirect_rh",i,"lh",j,sep="")
      indirect.eff[rh_i,lh_j]=if(!is.null(effects$indirect.effect.elem)){paste(indirect.eff.id[rh_i,lh_j],":=",paste0(effects$indirect.effect.elem,collapse ="+"),sep="")}else{""}

      direct.eff.id[rh_i,lh_j]=paste("direct_rh",i,"lh",j,sep="")
      direct.eff[rh_i,lh_j]=if(!is.null(effects$direct.effect.elem)){paste(direct.eff.id[rh_i,lh_j],":=",paste0(effects$direct.effect.elem,collapse ="+"),sep="")}else{""}
    }
  }
  all.effect=c(total.eff[],indirect.eff[],direct.eff[])
  all.effect=all.effect[!all.effect==""]
  effects=paste0(all.effect,collapse =";")

  model.with.effects=if(is.null(model)){NULL}else{paste(model,effects,sep=";")}

  mat[]=NA;mat[total.eff==""]=0;   mat[total.eff!=""]=1;   total.eff=mat
  mat[]=NA;mat[direct.eff==""]=0;  mat[direct.eff!=""]=1;  direct.eff=mat
  mat[]=NA;mat[indirect.eff==""]=0;mat[indirect.eff!=""]=1;indirect.eff=mat
  return(
    list(model=list(model.with.effects=model.with.effects, effects.defined=effects),
      total.effects=list(mat=t(total.eff),id=t(total.eff.id)),
      direct.effects=list(mat=t(direct.eff),id=t(direct.eff.id)),
      indirect.effects=list(mat=t(indirect.eff),id=t(indirect.eff.id))
    )
  )
}


