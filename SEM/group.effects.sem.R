group.effects.sem<-function(fit,rh=NULL,lh=NULL,model=NULL,
                            group.rhs=NULL,
                            group.lhs=NULL,
                            via=NULL,exclude.paths=NULL,exclude.pars=NULL){



  #effects.model(rh=c("D","I","PR"),lh=c("Multi.stibility","M_Wellbeing"),edges,model=MD)
  par.tab_model=extract.sem.model(fit)
  par.tab=par.tab_model$model.table
  edges=par.tab[par.tab$op=="~",]
  edges0=edges
  if(!is.null(exclude.paths)){
    judge_exclude.paths=edges$experession.nolable%in%gsub(pattern = " ",replacement = "",x = exclude.paths)
  }else{
    judge_exclude.paths=FALSE
  }
  if(!is.null(exclude.pars)){
    judge_exclude.paths=judge_exclude.paths|(edges$lhs%in%exclude.pars|edges$rhs%in%exclude.pars)
  }
  if(sum(judge_exclude.paths)){
    print(c(paste(sum(judge_exclude.paths),if(sum(judge_exclude.paths)==1){"path <"}else{"paths <"},
                  paste(edges$experession.nolable[judge_exclude.paths],collapse = " & "),
                  if(sum(judge_exclude.paths)==1){"> was "}else{"> were "},"not included in the calculation of the direct, indirect or total effects")))

    edges=edges[!judge_exclude.paths,]
  }

  if(is.null(rh)){
    rh=rownames(sortby(tabulate.char(edges$rhs),sort_vector = 1,decreasing = T))
  }
  if(is.null(lh)){
    lh=rownames(sortby(tabulate.char(edges$lhs),1))
  }

  j1=rh%in%lh
  if(sum(j1)>0&sum(j1)<length(rh)){rh=c(rh[!j1],rh[j1])}
  j2=lh%in%rh
  if(sum(j2)>0&sum(j2)<length(lh)){lh=c(lh[j2],lh[!j2])}

  rh.id=NULL
  for(i in 1:length(rh)){
    rh.id[rh[i]]=i
  }
  lh.id=NULL
  for(j in 1:length(lh)){
    lh.id[lh[j]]=j
  }


  if(!is.null(group.rhs)){
    if(is.character(group.rhs)){
      group.rhs=gsub(" ","",group.rhs)
      group.rhs=strsplit(x = group.rhs,split = "\\+")
    }
    temp=NULL
    sorunique=NULL
    equal.eff=NULL
    for(i in names(group.rhs)){
      lh.rh=sort(unique(group.rhs[[i]]))
      if(length(lh.rh)>1){
        lh.rh=lh.rh[!lh.rh %in% edges[edges$lhs%in%lh.rh & edges$rhs%in%lh.rh,]$lhs]
        if(length(lh.rh)>0){
          if(!paste0(lh.rh,collapse = "")%in%sorunique){
            temp[[i]]=lh.rh
          }else{
            equal.eff[[i]]=paste0(lh.rh,collapse = "")
          }

        }
      }else if(length(lh.rh)==1){
        if(!paste0(lh.rh,collapse = "")%in%sorunique){
          temp[[i]]=lh.rh
        }else{
          equal.eff[[i]]=paste0(lh.rh,collapse = "")
        }
      }
      if(!is.null(lh.rh)){sorunique=unique(c(sorunique,paste0(lh.rh,collapse = "")))}
      lh.rh=NULL
    }
    group.rhs=temp
  }else{
    group.rhs=as.list(rh);
    names(group.rhs)=rh
  }
  if(!is.null(equal.eff)){
    for(i in names(equal.eff)){
      log.i=as.matrix(lapply(group.rhs,FUN = paste0,collapse=""))
      equal.eff[[i]]=rownames(log.i)[equal.eff[[i]]==log.i]
    }
  }

  if(!is.null(group.lhs)){
    if(is.character(group.lhs)){
      group.lhs=gsub(" ","",group.lhs)
      group.lhs=strsplit(x = group.lhs,split = "\\+")
    }
    temp=NULL
    sorunique=NULL
    for(i in names(group.lhs)){
      lh.rh=sort(unique(group.lhs[[i]]))
      if(length(lh.rh)>0 & !paste0(lh.rh,collapse = "")%in%sorunique){
        temp[[i]]=lh.rh
      }
      if(!is.null(lh.rh)){sorunique=unique(c(sorunique,paste0(lh.rh,collapse = "")))}
      lh.rh=NULL
    }
    group.lhs=temp
  }else{
    group.lhs=as.list(lh);
    names(group.lhs)=lh
  }

  if(!is.null(via)){
    #if(length(names(via))<length(via)){
      #names(via)=c(names(via),paste("via",(length(via)-length(names(via))+1):length(via)))
    #}
    BY=NULL
    if(is.null(names(via))){names(via)=gsub(" ","",via)}
    for(via.i in names(via)){
      BY[[via.i]]=unlist(strsplit(gsub(" ","",via[via.i]),split = "\\+"))
    }
    by1=BY[[1]];names(by1)=by1;by_pars=NULL
    for(tempi in by1){
      by_pars[[tempi]]=edges$label[edges$lhs==tempi]
    }
  }else{
    by_pars=NULL
  }

  group.effect_mat=matrix(NA,length(group.lhs),length(group.rhs),dimnames = list(names(group.lhs),names(group.rhs)))
  group.total.effect_mat=group.direct.effect_mat=group.indirect.effect_mat=group.effect_mat
  group.effect_id=group.effect_mat
  group.effect_id[]=""
  group.total.effect_id=group.direct.effect_id=group.indirect.effect_id=group.effect_id
  if(!is.null(via)){
    indirect.by.mat=indirect.by.id=NULL
    for(i in names(by_pars)){
      indirect.by.mat[[i]]=group.effect_mat;indirect.by.id[[i]]=group.effect_id
    }
    indirect.by.mat[["interac"]]=group.effect_mat;indirect.by.id[["interac"]]=group.effect_id
    indirect.by.mat[["others"]]=group.effect_mat;indirect.by.id[["others"]]=group.effect_id

  }

  group.element=NULL
  for(i in names(group.rhs)){
    for (j in names(group.lhs)){
      id.rh=rh.id[group.rhs[[i]]]
      id.lh=lh.id[group.lhs[[j]]]
      temp=paste0("total_rh",id.rh,"lh",rep(id.lh,each=length(id.rh)))
      group.element[[paste0(i,"_VS_",j)]]=temp
    }
  }

  if(is.null(model)){
    model=par.tab_model$model
  }

  mat=matrix(NA,length(rh),length(lh),dimnames = list(rh,lh))
  total.eff=direct.eff=indirect.eff=NULL
  for(i in 1:length(rh)){
    rh_i=rh[i]
    for(j in 1:length(lh)){
      lh_j=lh[j]
      effects=SEM.total.effect.expression(rh=rh_i,lh=lh_j,edges)

      total.eff[[paste("total_rh",i,"lh",j,sep="")]]=effects$total.effect.elem
      indirect.eff[[paste("indirect_rh",i,"lh",j,sep="")]]=effects$indirect.effect.elem
      direct.eff[[paste("direct_rh",i,"lh",j,sep="")]]=effects$direct.effect.elem
    }
  }

  ###########################################################
  all.factor=unique(c(rh,lh))
  direct.path.ids=indirect.path.ids=NULL
  #group.pars=c(group.rhs,group.lhs[names(group.lhs)[!(names(group.lhs)%in%names(group.rhs))]])
  if(!is.null(via)){
    for(g.rhsi in names(group.element)){
      rhilhi=unlist(strsplit(g.rhsi,split = "_VS_"));names(rhilhi)=c("rh","lh")
      indirect.factors=BY[[1]]
      indirect.path.ids[[g.rhsi]]=edges$label[ (edges$lhs%in%indirect.factors | edges$rhs%in%indirect.factors)]
      direct.path.ids[[g.rhsi]]=edges$label[!(edges$lhs%in%indirect.factors | edges$rhs%in%indirect.factors)]
      "for(g.rhsj in names(group.element)[!names(group.element)==g.rhsi]){
        indirect.path.ids[[g.rhsj]]=sort(unique(c(indirect.path.ids[[g.rhsj]],direct.path.ids[[g.rhsi ]])))
       }"
    }
  }else{
    for(g.rhsi in names(group.element)){
      rhilhi=unlist(strsplit(g.rhsi,split = "_VS_"));names(rhilhi)=c("rh","lh")

      direct.factors=c(unlist(group.rhs[rhilhi["rh"]]),unlist(group.lhs[rhilhi["lh"]]))
      if(length(direct.factors)>0){if(length(direct.factors)>length(unique(direct.factors))){warning(paste(g.rhsi,": repeated pars in rhs and lhs"))}}
      direct.factors=unique(direct.factors)
      direct.path.ids[[g.rhsi]]=edges$label[ edges$lhs%in%direct.factors & edges$rhs%in%direct.factors ]

      indirect.factors=all.factor[!all.factor%in%direct.factors]
      indirect.path.ids[[g.rhsi]]=edges$label[ (edges$lhs%in%indirect.factors | edges$rhs%in%indirect.factors)]

    }
  }
  # g.rhsi="Eco.stability_VS_Eco.stability"

  alleffect.temp.total=sapply(names(total.eff),function(iii,...) paste0(total.eff[[iii]],collapse = "+"))

  def.total=def.direct=def.indirect=NULL;def.indirect.by=NULL

  for(g.rhsi in names(group.element)){
    rhilhi=unlist(strsplit(g.rhsi,split = "_VS_"));names(rhilhi)=c("rh","lh")
    direct.group.total=indirect.group.total=NULL
    direct.group.total1=indirect.group.total1=NULL
    indirect.group.by=NULL

    for( ii in group.element[[g.rhsi]] ){
      if(!is.null(total.eff[[ii]])){
        direct.ind=gsub("total_","direct_",ii)
        inderect.ind=gsub("total_","indirect_",ii)
        direct.ii=direct.eff[[direct.ind]]
        indirect.ii=indirect.eff[[inderect.ind]]
        total.ii=total.eff[[ii]]
        if(!all(sort(unique(total.ii))==sort(unique(c(direct.ii,indirect.ii))))){stop("codes error, please debug!!!")}

        with_indirect.path.pars=vector(length = length(indirect.ii))
        with_direct.path.pars=with_indirect.path.pars

        if(!is.null(indirect.ii)&&length(indirect.ii)>0&&!is.na(indirect.ii)){
          for(tempi in 1:length(indirect.ii)){
            indirect.ii.temp=unlist(strsplit(indirect.ii[tempi],split = "\\*"))
            with_indirect.path.pars[tempi]=any(indirect.ii.temp%in%indirect.path.ids[[g.rhsi]])
            with_direct.path.pars[tempi]=all(indirect.ii.temp%in%direct.path.ids[[g.rhsi]])
          }
          if(sum(with_indirect.path.pars)){
            indirect.group.total[[ii]]=indirect.ii[with_indirect.path.pars]
          }

          if(sum(!with_direct.path.pars)){
            indirect.group.total1[[ii]]=indirect.ii[!with_direct.path.pars]
          }
        }


        if(sum(!with_indirect.path.pars)){
          direct.group.total[[ii]]=c(direct.ii,indirect.ii[!with_indirect.path.pars])
        }else if(!is.null(direct.ii)&&length(direct.ii)>0&&!is.na(direct.ii)){
          direct.group.total[[ii]]=direct.ii
        }

        if(sum(with_direct.path.pars)){
          direct.group.total1[[ii]]=c(direct.ii,indirect.ii[with_direct.path.pars])
        }else if(!is.null(direct.ii)&&length(direct.ii)>0&&!is.na(direct.ii)){
          direct.group.total1[[ii]]=direct.ii
        }


        #if(!is.null(indirect.group.total1)){print(".....")}
      }
    }

    if(!all(unlist(indirect.group.total)==unlist(indirect.group.total1))|!all(unlist(direct.group.total)==unlist(direct.group.total1))){
      stop("line 179 in 'group.effects: unknown error !")
    }

    if(!is.null(direct.group.total)){
      direct.effect.temp=sapply(names(direct.group.total),function(i,...) paste0(direct.group.total[[i]],collapse = "+"))
      group.ele.i=direct.effect.temp[names(direct.effect.temp)%in%group.element[[g.rhsi]]]
      group.ele.i=group.ele.i[!is.na(group.ele.i)&group.ele.i!=""]
      if(!is.null(group.ele.i)&length(group.ele.i)>0){
        group.direct.effect_mat[rhilhi["lh"],rhilhi["rh"]]=1
        group.direct.effect_id[rhilhi["lh"],rhilhi["rh"]]=paste0("direct_",g.rhsi)
        def.direct=c(def.direct,paste(paste0("direct_",g.rhsi),paste0(group.ele.i,collapse = "+"),sep=":="))
      }

    }

    if(!is.null(indirect.group.total)){
      if(!is.null(via)){
        name.ind=names(indirect.group.total)
        name.ind.i=0
        for(tempi in indirect.group.total){
          name.ind.i=name.ind.i+1
          name=name.ind[name.ind.i]
          for(tempk in 1:length(tempi)){
            temp1=strsplit(tempi[tempk],"\\*")[[1]]
            judge_in=vector(length = length(by_pars));names(judge_in)=names(by_pars)
            for(tempj in names(by_pars)){
              judge_in[tempj]=any(temp1[1:(length(temp1)-1)]%in%by_pars[[tempj]])
            }

            if(sum(judge_in)==1){
              indirect.group.by[[names(judge_in)[judge_in]]][[name]]=
                c(indirect.group.by[[names(judge_in)[judge_in]]][[name]],tempi[tempk])
            }else if(sum(judge_in)>1){
              indirect.group.by[["interac"]][[name]]=c(indirect.group.by[["interac"]][[name]],tempi[tempk])
            }else{
              indirect.group.by[["others"]][[name]]=c(indirect.group.by[["others"]][[name]],tempi[tempk])
            }
          }
        }

        for(tempi in names(indirect.group.by)) {
          indirect.group.by.temp=sapply(names(indirect.group.by[[tempi]]),function(i,...) paste0(indirect.group.by[[tempi]][[i]],collapse = "+"))
          group.ele.i=indirect.group.by.temp[names(indirect.group.by.temp)%in%group.element[[g.rhsi]]]
          group.ele.i=group.ele.i[!is.na(group.ele.i)&group.ele.i!=""]
          if(!is.null(group.ele.i)&length(group.ele.i)>0){
            indirect.by.mat[[tempi]][rhilhi["lh"],rhilhi["rh"]]=1
            indirect.by.id[[tempi]][rhilhi["lh"],rhilhi["rh"]]=paste0("by.",tempi,"_",g.rhsi)
            def.indirect.by[[tempi]]=c(def.indirect.by[[tempi]],paste(paste0("by.",tempi,"_",g.rhsi),paste0(group.ele.i,collapse = "+"),sep=":="))
          }
        }
      }

      indirect.effect.temp=sapply(names(indirect.group.total),function(i,...) paste0(indirect.group.total[[i]],collapse = "+"))
      group.ele.i=indirect.effect.temp[names(indirect.effect.temp)%in%group.element[[g.rhsi]]]
      group.ele.i=group.ele.i[!is.na(group.ele.i)&group.ele.i!=""]
      if(!is.null(group.ele.i)&length(group.ele.i)>0){
        group.indirect.effect_mat[rhilhi["lh"],rhilhi["rh"]]=1
        group.indirect.effect_id[rhilhi["lh"],rhilhi["rh"]]=paste0("indirect_",g.rhsi)
        def.indirect=c(def.indirect,paste(paste0("indirect_",g.rhsi),paste0(group.ele.i,collapse = "+"),sep=":="))
      }
    }

    total.effect.total=total.eff[names(total.eff)%in%group.element[[g.rhsi]]]
    if(!is.null(total.effect.total)){
      total.effect.temp=sapply(names(total.effect.total),function(i,...) paste0(total.effect.total[[i]],collapse = "+"))
      group.ele.i=total.effect.temp[names(total.effect.temp)%in%group.element[[g.rhsi]]]
      group.ele.i=group.ele.i[!is.na(group.ele.i)&group.ele.i!=""]
      if(!is.null(group.ele.i)&length(group.ele.i)>0){
        group.total.effect_mat[rhilhi["lh"],rhilhi["rh"]]=1
        group.total.effect_id[rhilhi["lh"],rhilhi["rh"]]=paste0("total_",g.rhsi)
        def.total=c(def.total,paste(paste0("total_",g.rhsi),paste0(group.ele.i,collapse = "+"),sep=":="))
      }
    }

  }
  group.total.effect_mat[is.na(group.total.effect_mat)]=0
  group.direct.effect_mat[is.na(group.direct.effect_mat)]=0
  group.indirect.effect_mat[is.na(group.indirect.effect_mat)]=0

  effects.group.total=paste0(def.total,collapse =";")
  effects.group.direct=paste0(def.direct,collapse =";")
  effects.group.indirect=paste0(def.indirect,collapse =";")

  effects=paste0(c(effects.group.total,effects.group.direct,effects.group.indirect),collapse = ";")


  effects.group.indirect.by=NULL;indirect.effects.by=NULL
  if(!is.null(via)){
    for(i in names(def.indirect.by)){
      indirect.by.mat[[i]][is.na(indirect.by.mat[[i]])]=0
      if(!is.null(def.indirect.by[[i]])){effects.group.indirect.by[[i]]=paste0(def.indirect.by[[i]],collapse =";")}
      indirect.effects.by[[i]]$mat=indirect.by.mat[[i]]
      indirect.effects.by[[i]]$id=indirect.by.id[[i]]
    }
    indirect.by.mat=indirect.by.mat[names(def.indirect.by)]
    indirect.by.id=indirect.by.id[names(def.indirect.by)]



    effects.by=paste0(effects.group.indirect.by,collapse = ";")
  }else{
    effects.by=NULL
  }




  model.with.effects=if(is.null(model)){NULL}else{paste(model,effects,effects.by,sep=";")}




  return(
    list(model=list(model.with.effects=model.with.effects, effects.defined=effects),
         total.effects=list(mat=group.total.effect_mat,id=group.total.effect_id),
         direct.effects=list(mat=group.direct.effect_mat,id=group.direct.effect_id),
         indirect.effects=list(mat=group.indirect.effect_mat,id=group.indirect.effect_id),
         indirect.effects.by=indirect.effects.by,
         equal.eff=equal.eff)

  )


}

