grouped.sem.effects.stat=function(fit,data,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                                  group.rhs=NULL,
                                  group.lhs=NULL,
                                  via=NULL,
                                  rh=NULL,lh=NULL,standardized = TRUE,
                                  bootstrap=FALSE,reps=1000,
                                  return.boot=TRUE,missing="listwise",boot.depth=NULL,permutation=FALSE){
  
                                    
  MD.with.effects=group.effects.sem(fit,rh=NULL,lh=NULL,model=NULL,
                                    group.rhs=group.rhs,
                                    group.lhs=group.lhs,
                                    via=via)
  fit = lavaan(model = MD.with.effects$model$model.with.effects,
               data = data,missing=missing)
  Est <- parameterEstimates(fit, ci = FALSE, standardized = standardized,remove.system.eq = T,remove.nonfree=T,rsquare=T,fmi=F)
  effects.Est=Est[Est$op==":=",]
  rownames(effects.Est)=effects.Est$label
  
  
  effects=se.effects=p.effects=std_effects=se.std_effects=p.std_effects=judge.non.zero=NULL
  effects.summary=std.effects.summary=NULL
  if(!bootstrap){
    for(i in effect.type){
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
      return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                  std_effects=std_effects,
                  effects.summary=effects.summary,
                  std.effects.summary=std.effects.summary))
    }else{
      return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                  effects.summary=effects.summary))
    }
    # return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,std_effects=std_effects))
  }else{
    
    all.pars=unique(c(Est[!Est$op==":=",]$lhs,Est[!Est$op==":=",]$rhs))
    data.filter=data[,all.pars]
    effects_array=std_effects_array=NULL
    judge.non.zero=len=NULL
    
    for(i in effect.type){
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
    
    for(reps.i in 1:reps){
      boot.data<-if(is.null(boot.depth)){
        data.filter[sample.int(nrow(data.filter),replace = T),]
      }else{
        data.filter[sample.int(n = nrow(data.filter),size = boot.depth,replace = T),]
      }
      if(permutation){boot.data=apply(boot.data,MARGIN = 2,function(x) x[sample.int(length(x),replace = FALSE)])}
      fit = lavaan(model = MD.with.effects$model$model.with.effects,
                   data = boot.data)
      Est <- parameterEstimates(fit, ci = FALSE, standardized = standardized,remove.system.eq = T,remove.nonfree=T,rsquare=T,fmi=F)
      effects.Est=Est[Est$op==":=",]
      rownames(effects.Est)=effects.Est$label
      
      for(i in effect.type){
        est<-effects.Est[MD.with.effects[[i]]$id[judge.non.zero[[i]]],]
        
        effects_array[[i]][judge.non.zero[[i]] + (reps.i-1)*len[[i]]]=est$est
        if(standardized){
          std_effects_array[[i]][judge.non.zero[[i]]+(reps.i-1)*len[[i]]]=est$std.all
        }
      }
    }
    
    for(i in effect.type){
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
        return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                    std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                    effects.summary=effects.summary,
                    std.effects.summary=std.effects.summary,
                    std.effetc.boot=std_effects_array,
                    effects.boot=effects_array))
      }else{
        return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                    std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                    effects.summary=effects.summary,
                    std.effects.summary=std.effects.summary))
      }

    }else{
      if(return.boot){
        return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                    effects.summary=effects.summary,
                    effects.boot=effects_array))
      }else{
        return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                    effects.summary=effects.summary))
      }
 
    }
    
  }
  
  
}

