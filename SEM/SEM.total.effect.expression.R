SEM.total.effect.expression<-function(rh=NULL,lh=NULL,edges){
  #total.effect.expression(rh="D",lh="Multi.stibility",edges)
  fun1<-function(x,edges){
    "x=sortby(x,-(1:length(x)))";
    if(length(x)<2&length(unique(x))<length(x)){
      return(NULL)
    }else if(length(x)==2){
      return(edges[paste0(x[2:1],collapse = "~"),]$label)
    }else if(length(x)>2){
      temp=vector(length = length(x)-1)
      for(k in 2:length(x)){
        temp[k-1]=paste0(x[k:(k-1)],collapse ="~")
      }
      return(paste0(edges[temp,]$label,collapse="*"))
    }
  }

  if(sum(edges$lhs==lh)==0|sum(edges$rhs==rh)==0){
    finish=NULL
  }else if(sum(edges$lhs==lh)==1|sum(edges$rhs==rh)==1){
    if(!paste0(lh,"~",rh)%in%edges$experession.nolable){
      finish=NULL
    }else{
      finish=paste0(rh,"~",lh)
    }
  }else{
    rh.temp=rh
    finish=NULL
    unfinish.member=rh
    n=0
    while (length(rh.temp)>0&n<(length(unique(edges$lhs)))){
      n=n+1
      unfinish.rh=NULL
      unfinish=NULL
      for(rhi in unique(rh.temp)){
        times.rhi=rh.temp==rhi

        lh.temp=edges$lhs[edges$rhs==rhi]
        if(length(lh.temp)>0){
          lh.temp=rep(lh.temp,each=sum(times.rhi))
          find.lh=lh.temp==lh
          #find.lh=rep(find.lh,times=sum(times.rhi))
          temp1=paste(unfinish.member[times.rhi],lh.temp,sep = "~")
          rh.new=if(sum(!find.lh)){lh.temp[!find.lh]}else{NULL}
          if(sum(find.lh)){
            finish=c(finish,temp1[find.lh]);
            temp1=temp1[!find.lh]
          }
          unfinish=c(unfinish,temp1)
          unfinish.rh=c(unfinish.rh,rh.new)
        }else{
          #print(rhi)
          next
        }

      }
      rh.temp=unfinish.rh
      unfinish.member=unfinish
    }

  }



  if(!is.null(finish)){
    finish.split=strsplit(x =finish,split = "~" )

    total.effect.elem=sapply(1:length(finish.split),function(i,...) fun1(finish.split[[i]],edges))

    direct_edge=paste(lh,rh,sep="~")==rownames(edges)
    direct.effect.elem=if(sum(direct_edge)){edges[direct_edge,]$label}else{NULL}
    indirect.effect.elem=
      if(is.null(direct.effect.elem)){
        total.effect.elem
      }else{
        if(length(total.effect.elem)>1){total.effect.elem[!total.effect.elem==direct.effect.elem]}else{NULL}
      }

  }else{
    total.effect.elem=direct.effect.elem=indirect.effect.elem=NULL
  }


  return(list(total.effect.elem=total.effect.elem,
              direct.effect.elem=direct.effect.elem,
              indirect.effect.elem=indirect.effect.elem
  ))
}



