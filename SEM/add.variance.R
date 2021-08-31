add.variance=function(model=NULL,variance=NULL,covariance=NULL){
  if(!is.null(model)){
    temp=strsplit(model,split = "[\n;]");
    md=NULL
    for(i in 1:length(temp)){
      md=c(md,temp[[i]])
    }
    md=gsub("[ ]*","",md)
    md=md[!grepl("^#",md)&!md==""]
    all.pars=unique(unlist(strsplit(md,"[~=:\\+]")));all.pars=all.pars[!all.pars==""]
    model.list=list(variance=NULL,covariance=NULL,effect=NULL,lv=NULL,def.v=NULL)
    if(sum(grepl("=~",md))){model.list$lv=md[grepl("=~",md)];md=md[!grepl("=~",md)]}

    if(sum(grepl(":=",md))){model.list$def.v=md[grepl(":=",md)];  md=md[!grepl(":=",md)]}

    if(sum(grepl("~~",md))){
      temp=md[grepl("~~",md)];
      md=md[!grepl("~~",md)]
      temp1=strsplit(temp,"~~")
      for(i in 1:length(temp1)){
        tempj=strsplit(temp1[[i]],"\\+")
        for(j1 in tempj[[1]]){
          for(j2 in tempj[[2]]){
            if(j1==j2){
              model.list$variance=c(model.list$variance,paste0(j1,"~~",j2))
            }else{
              model.list$covariance=c(model.list$covariance,paste0(j1,"~~",j2))
            }
          }
        }
      }
    }
    md=md[grepl("~",md)]
    temp1=strsplit(md,"~")
    if(length(temp1)>0){
      for(i in 1:length(temp1)){
        tempj=strsplit(temp1[[i]],"\\+")
        for(j1 in tempj[[1]]){
          for(j2 in tempj[[2]]){
            model.list$effect=c(model.list$effect,paste0(j1,"~",j2))
          }
        }
      }
    }

  }
  if(is.null(variance)){variance=all.pars}
  variance=unique(c(variance,all.pars))
  model.list$variance=unique(c(model.list$variance,paste0(variance,"~~",variance)))
  tab=t(as.data.frame(strsplit(x = unlist(model.list),split = "~[~]*")))
  if(length(covariance)>0){
    covariance=unique(covariance)
    for(j1 in 1:length(covariance)){
      for (j2 in j1:length(covariance)){
        if(j1==j2){
          model.list$variance=c(model.list$variance,paste0(covariance[j1],"~~",covariance[j2]))
        }else{
          model.list$covariance=c(model.list$covariance,paste0(covariance[j1],"~~",covariance[j2]))
        }
      }
    }
    model.list$covariance=unique(sort(model.list$covariance))
    model.list$covariance=model.list$covariance[!gsub("~","",model.list$covariance)%in%c(paste0(tab[,1],tab[,2]),paste0(tab[,2],tab[,1]))]
    model.list$variance=unique(sort(model.list$variance))
  }

  for(i in 1:length(model.list)){
    if(!is.null(model.list[[i]])){model.list[[i]]=sort(unique(model.list[[i]]))}
  }
  model=paste0(unique(unlist(model.list)),collapse = ";")
  return(model)
}
