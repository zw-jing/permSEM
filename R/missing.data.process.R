missing.data.process <- function(x,method=c("rnorm","resample.self"),replace=TRUE){
  #par(mfrow=c(2,2))
  #x=c(rnorm(10000,mean = -5,sd = 1),rnorm(10000,mean =5,sd = 1))
  #hist(x,breaks = 100)
  #x[sample.int(20000,10000)]=NA
  #hist(x[!is.na(x)],breaks = 100)
  #hist(missing.data.process(x,method = "rnorm"),breaks = 100)
  #hist(missing.data.process(x,method = "re"),breaks = 100)
  x[is.na(x)]=NA
  na.logic=is.na(x)
  if(sum(na.logic)>0&sum(!na.logic)){
    if(method[1] == "rnorm"){
      rand.norm.set=rnorm(sum(na.logic)*10,mean(x,na.rm=T),sd(x,na.rm=T))
      j=rand.norm.set<=max(x[!na.logic])&rand.norm.set>=min(x[!na.logic])
      x[na.logic]=sample(rand.norm.set[j],sum(na.logic))
    }else{
      rand.norm.set=spline(x[!na.logic],n = sum(!na.logic)*(ceiling(sum(na.logic)/sum(!na.logic))+2))
      j=rand.norm.set$y<=max(x[!na.logic])&rand.norm.set$y>=min(x[!na.logic])
      x[na.logic]=sample(rand.norm.set$y[j],sum(na.logic),replace = replace)
    }
  }
  return(x)
}


