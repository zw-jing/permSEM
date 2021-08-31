grouped.sem.effects.stat.permute_strata<-function(fit,data,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                                           group.rhs=NULL,
                                           group.lhs=NULL,
                                           via=NULL,
                                           rh=NULL,lh=NULL,standardized = TRUE,
                                           bootstrap=FALSE,reps=1000,
                                           return.boot=TRUE,missing="listwise",boot.depth=NULL,permutation=FALSE,
                                           ind.in.dist=NULL,dist.list=NULL,boot.perm.simultaneous = FALSE,
                                           strata=NULL,exclude.paths=NULL,exclude.pars=NULL){

  if(!is.null(ind.in.dist)&length(ind.in.dist)!=nrow(data)){stop("ind.in.dist input wrong")}
  MD.with.effects=group.effects.sem(fit,rh=NULL,lh=NULL,model=NULL,
                                    group.rhs=group.rhs,
                                    group.lhs=group.lhs,
                                    via=via,exclude.paths=exclude.paths,exclude.pars=exclude.pars)
  if(!is.null(MD.with.effects$indirect.effects.by)){
    indirect.by_ele=paste0("indirect.by_",names(MD.with.effects$indirect.effects.by))
    for(i in names(MD.with.effects$indirect.effects.by)){
      MD.with.effects[[paste0("indirect.by_",i)]]=MD.with.effects$indirect.effects.by[[i]]
    }
    eff.type=c(effect.type,indirect.by_ele)
  }
  fit = lavaan(model = MD.with.effects$model$model.with.effects,
               data = data,missing=missing)
  Est <- parameterEstimates(fit, ci = FALSE, standardized = standardized,remove.system.eq = T,remove.nonfree=T,rsquare=T,fmi=F)
  effects.Est=Est[Est$op==":=",]
  rownames(effects.Est)=effects.Est$label

  if(permutation){boot.depth=NULL}
  if(!bootstrap){boot.perm.simultaneous=FALSE}

  effects=se.effects=p.effects=std_effects=se.std_effects=p.std_effects=judge.non.zero=NULL
  effects.summary=std.effects.summary=NULL
  if(!bootstrap){
    for(i in eff.type){
      temp.mat=MD.with.effects[[i]]$mat
      judge.non.zero[[i]]=which(temp.mat>0)
      temp.mat[]=NA
      effects[[i]]=se.effects[[i]]=p.effects[[i]]=temp.mat


      est=effects.Est[MD.with.effects[[i]]$id[judge.non.zero[[i]]],]

      effects[[i]][judge.non.zero[[i]]]=est$est
      se.effects[[i]][judge.non.zero[[i]]]=est$se
      p.effects[[i]][judge.non.zero[[i]]]=est$pvalue
      effects[[i]][is.na(effects[[i]])]=0

      effects.summary[[i]]=temp.mat
      summary.temp=paste(char.align(effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
      effects.summary[[i]][]=summary.temp[]


      if(standardized){
        std_effects[[i]]=se.std_effects[[i]]=p.std_effects[[i]]=temp.mat
        std_effects[[i]][judge.non.zero[[i]]]=est$std.all
        std_effects[[i]][is.na(std_effects[[i]])]=0

        std.effects.summary[[i]]=temp.mat
        summary.temp=paste(char.align(std_effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
        std.effects.summary[[i]][]=summary.temp[]
      }
    }


    if(standardized){
      return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                            std_effects=std_effects,
                            effects.summary=effects.summary,
                            std.effects.summary=std.effects.summary)
    }else{
      return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                            effects.summary=effects.summary)
    }
    # return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,std_effects=std_effects))
  }else{

    all.pars=unique(c(Est[!Est$op==":=",]$lhs,Est[!Est$op==":=",]$rhs))
    data.filter=data[,all.pars]
    effects_array=std_effects_array=NULL
    judge.non.zero=len=NULL

    for(i in eff.type){
      temp.mat=MD.with.effects[[i]]$mat;temp.mat[temp.mat==0]=NA
      effects_array[[i]]=array(NA,dim = c(nrow(temp.mat),ncol(temp.mat),reps))
      judge.non.zero[[i]]=which(MD.with.effects[[i]]$mat>0)
      len[[i]]=length(MD.with.effects[[i]]$mat)
      effects[[i]]=se.effects[[i]]=p.effects[[i]]=temp.mat

      effects.summary[[i]]=temp.mat;effects.summary[[i]][]=NA

      if(standardized){
        std_effects_array[[i]]=effects_array[[i]]
        std_effects[[i]]=se.std_effects[[i]]=p.std_effects[[i]]=temp.mat
        std.effects.summary[[i]]=effects.summary[[i]]
      }
    }
    re.samp.n=0
    for(reps.i in 1:reps){

      fit=NULL
      boot.data<-if(is.null(boot.depth)){
        data.filter[sample.int(nrow(data.filter),replace = T),]
      }else{
        data.filter[sample.int(n = nrow(data.filter),size = boot.depth,replace = T),]
      }
      #if(permutation){boot.data=apply(boot.data,MARGIN = 2,function(x) x[sample.int(length(x),replace = FALSE)])}
      fit = lavaan(model = MD.with.effects$model$model.with.effects,
                   data = boot.data)
      while(is.null(fit)){
        re.samp.n=re.samp.n+1
        boot.data<-if(is.null(boot.depth)){
          data.filter[sample.int(nrow(data.filter),replace = T),]
        }else{
          data.filter[sample.int(n = nrow(data.filter),size = boot.depth,replace = T),]
        }
        #if(permutation){boot.data=apply(boot.data,MARGIN = 2,function(x) x[sample.int(length(x),replace = FALSE)])}
        fit = lavaan(model = MD.with.effects$model$model.with.effects,
                     data = boot.data)
      }

      Est <- parameterEstimates(fit, ci = FALSE, standardized = standardized,remove.system.eq = T,remove.nonfree=T,rsquare=T,fmi=F)
      effects.Est=Est[Est$op==":=",]
      rownames(effects.Est)=effects.Est$label

      for(i in eff.type){
        est<-effects.Est[MD.with.effects[[i]]$id[judge.non.zero[[i]]],]

        effects_array[[i]][judge.non.zero[[i]] + (reps.i-1)*len[[i]]]=est$est
        if(standardized){
          std_effects_array[[i]][judge.non.zero[[i]]+(reps.i-1)*len[[i]]]=est$std.all
        }
      }
    }
    if(re.samp.n>0){warning(paste("Bootstraping waring:",re.samp.n,"moedles without solutions was recalculated"))}
    for(i in eff.type){
      effects[[i]][]=apply(X = effects_array[[i]],MARGIN = c(1,2),FUN = mean)
      se.effects[[i]][]=apply(X = effects_array[[i]],MARGIN = c(1,2),FUN = sd)
      p.effects[[i]][judge.non.zero[[i]]]=pt(-abs(effects[[i]][judge.non.zero[[i]]])/(se.effects[[i]][judge.non.zero[[i]]]+10^-26),df = reps-1)
      effects[[i]][is.na(effects[[i]])]=0

      summary.temp=paste(char.align(effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
      effects.summary[[i]][]=summary.temp[]

      if(standardized){
        std_effects[[i]][]=apply(X = std_effects_array[[i]],MARGIN = c(1,2),FUN = mean)
        se.std_effects[[i]][]=apply(X = std_effects_array[[i]],MARGIN = c(1,2),FUN = sd)
        p.std_effects[[i]][judge.non.zero[[i]]]=pt(-abs(std_effects[[i]][judge.non.zero[[i]]])/(se.std_effects[[i]][judge.non.zero[[i]]]+10^-26),df = reps-1)
        std_effects[[i]][is.na(std_effects[[i]])]=0

        summary.temp=paste(char.align(std_effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.std_effects[[i]],align = "left",digits = 2))
        std.effects.summary[[i]][]=summary.temp[]
      }
    }

    if(standardized){
      if(return.boot){
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                              effects.summary=effects.summary,
                              std.effects.summary=std.effects.summary,
                              std.effetc.boot=std_effects_array,
                              effects.boot=effects_array)
      }else{
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                              effects.summary=effects.summary,
                              std.effects.summary=std.effects.summary)
      }

    }else{
      if(return.boot){
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              effects.summary=effects.summary,
                              effects.boot=effects_array)
      }else{
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              effects.summary=effects.summary)
      }

    }

  }

  if(permutation){
    all.pars=unique(c(Est[!Est$op==":=",]$lhs,Est[!Est$op==":=",]$rhs))
    perm.data=data[,all.pars]
    independ.pars=all.pars[!all.pars%in%Est[Est$op=="~",]$lhs]
    independ.df=data[,independ.pars]
    depend.pars=all.pars[!all.pars%in%independ.pars]
    depend.df=data[,depend.pars]
    perm.list=dist.list[depend.pars]


    temp.distmat=perm.list[[1]];temp.distmat[]<-0
    temp.distmat[ind.in.dist]=1
    subset.perm=sort(unique(c(which(apply(temp.distmat,1,sum)>0),which(apply(temp.distmat,2,sum)>0))))
    ind.in.dist.new=which(temp.distmat[subset.perm,subset.perm]>0)
    if(!all(temp.distmat[subset.perm,subset.perm][ind.in.dist.new]==temp.distmat[ind.in.dist])){
      stop("ind.in.dist parameter wrong!")
    }
    if(nrow(temp.distmat)>length(subset.perm)){
      ind.in.dist=ind.in.dist.new
      for(iii in 1:length(perm.list)){
        perm.list[[iii]]=perm.list[[iii]][subset.perm,subset.perm]
      }
    }

    effects_array=std_effects_array=NULL
    judge.non.zero=len=NULL

    for(i in eff.type){
      temp.mat=MD.with.effects[[i]]$mat;temp.mat[temp.mat==0]=NA
      effects_array[[i]]=array(NA,dim = c(nrow(temp.mat),ncol(temp.mat),reps))
      judge.non.zero[[i]]=which(MD.with.effects[[i]]$mat>0)
      len[[i]]=length(MD.with.effects[[i]]$mat)
      effects[[i]]=se.effects[[i]]=p.effects[[i]]=temp.mat

      effects.summary[[i]]=temp.mat;effects.summary[[i]][]=NA

      if(standardized){
        std_effects_array[[i]]=effects_array[[i]]
        std_effects[[i]]=se.std_effects[[i]]=p.std_effects[[i]]=temp.mat
        std.effects.summary[[i]]=effects.summary[[i]]
      }
    }

    #perm.mat=sapply(1:reps,function(iii,...) {sample.int(nrow(perm.list[[1]]))})
    perm.mat=t(vegan:::getPermuteMatrix(N = nrow(perm.list[[1]]),perm = min(reps*10,100000),strata = strata))
    perm.vec=vector(length = length(perm.list));names(perm.vec)=names(perm.list)
    len_perm=ncol(perm.mat)
    re.samp.n=0
    for(reps.i in 1:reps){
      print(reps.i)
      perm.vec[]<-sample.int(len_perm,length(perm.list))
      for(k in names(perm.list)){
        perm.data[[k]]=perm.list[[k]][perm.mat[,perm.vec[k]],perm.mat[,perm.vec[k]]][ind.in.dist]
      }
      if(standardized){
        for(tempi in depend.pars){
          perm.data[,tempi]<-(perm.data[,tempi]-mean(perm.data[,tempi],na.rm=T))/sd(perm.data[,tempi],na.rm=T)
        }
      }
      fit=NULL
      fit = lavaan(model = MD.with.effects$model$model.with.effects,
                   data = perm.data)

      while(is.null(fit)){
        re.samp.n=re.samp.n+1
        perm.vec=sample.int(nrow(perm.list[[1]]))
        for(k in names(perm.list)){
          perm.data[[k]]=perm.list[[k]][perm.vec,perm.vec][ind.in.dist]
        }
        if(standardized){
          for(tempi in depend.pars){
            perm.data[,tempi]<-(perm.data[,tempi]-mean(perm.data[,tempi],na.rm=T))/sd(perm.data[,tempi],na.rm=T)
          }
        }
        fit = lavaan(model = MD.with.effects$model$model.with.effects,
                     data = perm.data)
      }

      Est <- parameterEstimates(fit, ci = FALSE, standardized = standardized,remove.system.eq = T,remove.nonfree=T,rsquare=T,fmi=F)
      effects.Est=Est[Est$op==":=",]
      rownames(effects.Est)=effects.Est$label

      for(i in eff.type){
        est<-effects.Est[MD.with.effects[[i]]$id[judge.non.zero[[i]]],]

        effects_array[[i]][judge.non.zero[[i]] + (reps.i-1)*len[[i]]]=est$est
        if(standardized){
          std_effects_array[[i]][judge.non.zero[[i]]+(reps.i-1)*len[[i]]]=est$std.all
        }
      }
    }
    if(re.samp.n>0){warning(paste("Permutation waring:",re.samp.n,"moedles without solutions was recalculated"))}
    for(i in eff.type){
      if(boot.perm.simultaneous) {
        effects_array[[i]]=-effects_array[[i]]+return.result$effects.boot[[i]]
      }else{
        effects_array[[i]]=-effects_array[[i]]+rep(return.result$effects[[i]],dim(effects_array[[i]])[3])
      }
      effects[[i]][]=apply(X = effects_array[[i]],MARGIN = c(1,2),FUN = mean)
      se.effects[[i]][]=apply(X = effects_array[[i]],MARGIN = c(1,2),FUN = sd)
      p.effects[[i]][judge.non.zero[[i]]]=pt(-abs(effects[[i]][judge.non.zero[[i]]])/(se.effects[[i]][judge.non.zero[[i]]]+10^-26),df = reps-1)
      effects[[i]][is.na(effects[[i]])]=0

      summary.temp=paste(char.align(effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
      effects.summary[[i]][]=summary.temp[]

      if(standardized){
        if(boot.perm.simultaneous) {
          std_effects_array[[i]] = -std_effects_array[[i]] + return.result$std.effetc.boot[[i]]
        }else{
          std_effects_array[[i]] = -std_effects_array[[i]]+rep(return.result$effects[[i]],dim(std_effects_array[[i]])[3])
        }
        std_effects[[i]][]=apply(X = std_effects_array[[i]],MARGIN = c(1,2),FUN = mean)
        se.std_effects[[i]][]=apply(X = std_effects_array[[i]],MARGIN = c(1,2),FUN = sd)
        p.std_effects[[i]][judge.non.zero[[i]]]=pt(-abs(std_effects[[i]][judge.non.zero[[i]]])/(se.std_effects[[i]][judge.non.zero[[i]]]+10^-26),df = reps-1)
        std_effects[[i]][is.na(std_effects[[i]])]=0

        summary.temp=paste(char.align(std_effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.std_effects[[i]],align = "left",digits = 2))
        std.effects.summary[[i]][]=summary.temp[]
      }
    }

    if(standardized){
      if(return.boot){
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                              effects.summary=effects.summary,
                              std.effects.summary=std.effects.summary,
                              std.effetc.boot=std_effects_array,
                              effects.boot=effects_array)
      }else{
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                              effects.summary=effects.summary,
                              std.effects.summary=std.effects.summary)
      }

    }else{
      if(return.boot){
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              effects.summary=effects.summary,
                              effects.boot=effects_array)
      }else{
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              effects.summary=effects.summary)
      }

    }
  }

  rown=rownames(return.result$effects$total.effects)
  coln=colnames(return.result$effects$total.effects)
  for(i in c("std.effetc.boot","effects.boot")){
    if(i %in% names(return.result)){
      for(j in 1:length(return.result[[i]])){
        rownames(return.result[[i]][[j]])=rown
        colnames(return.result[[i]][[j]])=coln
      }
    }
  }



  #########################################################################################################
  if(!is.null(MD.with.effects$equal.eff)){
    equal.eff=MD.with.effects$equal.eff
    warning(paste(c("The effects of <",names(equal.eff),"> is equal to the effects of <",unlist(equal.eff),">"),collapse = " "))

    temp.array=array(NA,dim = c(length(rown),length(coln)+length(equal.eff),reps),dimnames = list(rown,c(coln,names(equal.eff)),1:reps))
    temp.mat=temp.array[,,1]
    cycle=data=NULL
    loop.depth=1
    for(cycle_1 in 1:length(return.result)){
      data[[loop.depth]]=return.result[[cycle_1]]


      if(!(is.array(data[[loop.depth]])|is.matrix(data[[loop.depth]])|is.data.frame(data[[loop.depth]]))){
        loop.depth=loop.depth+1
        for(cycle_2 in 1:length(data[[loop.depth-1]])){
          data[[loop.depth]]=data[[loop.depth-1]][[cycle_2]]
          if(!(is.array(data[[loop.depth]])|is.matrix(data[[loop.depth]])|is.data.frame(data[[loop.depth]]))){
            loop.depth=loop.depth+1
            warning("debug ... loop.depth low")
          }else{
            if(length(dim(return.result[[cycle_1]][[cycle_2]]))==3){
              temp.array[,coln,]=return.result[[cycle_1]][[cycle_2]]
              temp.array[,names(equal.eff),]=return.result[[cycle_1]][[cycle_2]][,unlist(equal.eff),]
              return.result[[cycle_1]][[cycle_2]]=temp.array
            }else if(length(dim(return.result[[cycle_1]][[cycle_2]]))==2){
              temp.mat[,coln]=return.result[[cycle_1]][[cycle_2]]
              temp.mat[,names(equal.eff)]=return.result[[cycle_1]][[cycle_2]][,unlist(equal.eff)]
              return.result[[cycle_1]][[cycle_2]]=temp.mat
            }

          }

        }
        loop.depth=loop.depth-1

      }else{
        if(length(dim(return.result[[cycle_1]]))==3){
          temp.array[,coln,]=return.result[[cycle_1]]
          temp.array[,names(equal.eff),]=return.result[[cycle_1]][,unlist(equal.eff),]
          return.result[[cycle_1]]=temp.array
        }else if(length(dim(return.result[[cycle_1]]))==2){
          temp.mat[,coln]=return.result[[cycle_1]]
          temp.mat[,names(equal.eff)]=return.result[[cycle_1]][,unlist(equal.eff)]
          return.result[[cycle_1]]=temp.mat
        }
      }
    }
  }

  #gc()
  return(return.result)

}

