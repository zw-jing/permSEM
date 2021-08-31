"CI.fun" <- function(data,confidance=0.05){
  # CI.fun(rnorm(1000))
  data <- as.matrix(data)
  fun.temp <- function(d, confidance=0.05){d=d[!is.na(d)];d=sort(d);len=length(d);n=confidance/2*len;w=(ceiling(n)-n);return(c(lower=(1-w)*d[ceiling(n)]+(w)*d[floor(n)],upper=w*d[ceiling(len-n)]+(1-w)*d[floor(len-n)]))}
  ci <- if(ncol(as.matrix(data)>1)){apply(data,2,fun.temp, confidance=confidance)}else{fun.temp(data, confidance=confidance)}
  return(ci)
}
"pool.CI.Fun" <- function(data.mat,MARGIN=2, permutations=1000, confidance=0.05,na.resample=TRUE, na.rm = TRUE){
  # pool.CI.Fun(rnorm(1000))
  pool.CI <- function(perm.data, d, N, confidance = 0.05, ...){
    perm.data[]=d[perm.data]
    if(na.resample){isna=is.na(perm.data);if(sum(isna)>0){perm.data[isna]=sample(d[!is.na(d)],sum(isna))}}
    return(CI=apply(t(CI.fun(perm.data)),2,mean,na.rm = na.rm))
  }

  data=as.matrix(data.mat)
  if(MARGIN == 1){data = t(data)}
  N=nrow(data)
  perm_mat=matrix(sample.int(N,size = N*permutations,replace = T),nrow=N,ncol = permutations)

  pooled.CI=if(ncol(data)==1){pool.CI(perm_mat, as.matrix(data[,1]), N=N, confidance = confidance)}else{sapply(1:ncol(data),function(i, ...) pool.CI(perm_mat, data[,i], N=N, confidance = confidance))}
  pool.CI=if(ncol(data)==1){as.matrix(as.data.frame(pooled.CI))}else{as.matrix(as.data.frame(pooled.CI))}
  colnames(pool.CI)=if(!is.null(colnames(data))){colnames(data)}else(1:ncol(data))
  return(pool.CI)
}
"pool.CI.ofMEAN.Fun" <- function(data.mat,MARGIN=2, permutations=1000, confidance=0.05,na.resample=TRUE, na.rm = TRUE){
  # pool.CI.ofMEAN.Fun(rnorm(1000))
  data=as.matrix(data.mat)
  if(MARGIN == 1){data = t(data)}
  N=nrow(data)
  perm_mat <- matrix(sample.int(N,size = N*permutations,replace = T),nrow=N,ncol = permutations)
  pool.ciofMean <- function(d, perm.data = NULL,permutations = 1000,confidance = 0.05){
    N <- length(d)
    if(is.null(perm.data)){
      perm.data <- matrix(sample.int(N,size = N*permutations,replace = TRUE),nrow=N,ncol = permutations)
    }
    perm.data[]=d[perm.data]
    if(na.resample){isna=is.na(perm.data);if(sum(isna)>0){perm.data[isna]=sample(d[!is.na(d)],sum(isna))}}
    return(CI.fun(mean.fun(perm.data),confidance = confidance))
  }
  pooled.CIofmean=if(ncol(data)==1){pool.ciofMean(as.matrix(data[,1]), perm_mat, permutations, confidance)}else{sapply(1:ncol(data),function(i, ...) pool.ciofMean(data[,i], perm_mat, permutations, confidance))}

  colnames(pooled.CIofmean)=colnames(data)
  rownames(pooled.CIofmean)=c("lower","upper")
  return(pooled.CIofmean)
}
"sd.fun"  <- function(data,na.rm=TRUE){if(ncol(as.matrix(data)>1)){apply(data,2,sd,na.rm=na.rm)}else{sd(data,na.rm=na.rm)}}
"mean.fun" <- function(data,na.rm=TRUE){if(ncol(as.matrix(data)>1)){apply(data,2,mean,na.rm=na.rm)}else{mean(data,na.rm=na.rm)}}
"pool.Sd.Fun" <- function(data.mat, MARGIN = 2, permutations=1000, na.resample=TRUE, na.rm=TRUE){
  # pool.Sd.Fun(rnorm(1000))
  data=as.matrix(data.mat)
  if(MARGIN == 1){data = t(data)}
  N=nrow(data)
  perm_mat <- matrix(sample.int(N,size = N*permutations,replace = T),nrow=N,ncol = permutations)
  pool.Sd <- function(d, perm.data = NULL,permutations = 1000){
    N <- length(d)
    if(is.null(perm.data)){
      perm.data <- matrix(sample.int(N,size = N*permutations,replace = TRUE),nrow=N,ncol = permutations)
    }
    perm.data[]=d[perm.data]
    if(na.resample){isna=is.na(perm.data);if(sum(isna)>0){perm.data[isna]=sample(d[!is.na(d)],sum(isna))}}
    return(mean(sd.fun(perm.data)))
  }
  pool.SD=if(ncol(data)==1){pool.Sd(as.matrix(data[,1]),perm_mat,permutations)}else{sapply(1:ncol(data),function(i, ...) pool.Sd(data[,i],perm_mat,permutations))}
  names(pool.SD)=colnames(data)
  return(pool.SD)
}
"pool.Mean.Fun" <- function(data.mat, MARGIN = 2, permutations=1000, na.resample=TRUE, na.rm=TRUE){
  #pool.Mean.Fun(rnorm(1000))
  data=as.matrix(data.mat)
  if(MARGIN == 1){data = t(data)}
  N=nrow(data)
  perm_mat <- matrix(sample.int(N,size = N*permutations,replace = T),nrow=N,ncol = permutations)
  pool.Mean <- function(d, perm.data = NULL,permutations = 1000){
    N <- length(d)
    if(is.null(perm.data)){
      perm.data <- matrix(sample.int(N,size = N*permutations,replace = TRUE),nrow=N,ncol = permutations)
    }
    perm.data[]=d[perm.data]
    if(na.resample){isna=is.na(perm.data);if(sum(isna)>0){perm.data[isna]=sample(d[!is.na(d)],sum(isna))}}
    return(mean(mean.fun(perm.data)))
  }
  pool.MEAN=if(ncol(data)==1){pool.Mean(as.matrix(data[,1]),perm_mat,permutations)}else{sapply(1:ncol(data),function(i, ...) pool.Mean(data[,i],perm_mat,permutations))}
  names(pool.MEAN)=colnames(data)
  return(pool.MEAN)
}
"pool.SE.Fun" <- function(data.mat, MARGIN = 2, permutations=1000, na.resample=TRUE, na.rm=TRUE){
  #pool.SE.Fun(rnorm(1000))
  data=as.matrix(data.mat)
  if(MARGIN == 1){data = t(data)}
  N=nrow(data)
  perm_mat <- matrix(sample.int(N,size = N*permutations,replace = T),nrow=N,ncol = permutations)
  pool.Sd.of.Mean <- function(d, perm.data = NULL,permutations = 1000){
    N <- length(d)
    if(is.null(perm.data)){
      perm.data <- matrix(sample.int(N,size = N*permutations,replace = TRUE),nrow=N,ncol = permutations)
    }
    perm.data[]=d[perm.data]
    if(na.resample){isna=is.na(perm.data);if(sum(isna)>0){perm.data[isna]=sample(d[!is.na(d)],sum(isna))}}
    return(sd(mean.fun(perm.data)))
  }
  pool.SE=if(ncol(data)==1){pool.Sd.of.Mean(as.matrix(data[,1]),perm_mat,permutations)}else{sapply(1:ncol(data),function(i, ...) pool.Sd.of.Mean(data[,i],perm_mat,permutations))}
  names(pool.SE)=colnames(data)
  return(pool.SE)
}
"pool.mean_sd_se" <- function(data.mat, MARGIN = 2, permutations=1000, na.resample=TRUE, na.rm=TRUE){
  #pool.mean_sd_se(rnorm(1000))

  data=as.matrix(data.mat)
  if(MARGIN == 1){data = t(data)}
  N=nrow(data)
  perm_mat <- matrix(sample.int(N,size = N*permutations,replace = T),nrow=N,ncol = permutations)
  fun.mix <- function(d, perm.data = NULL,permutations = 1000){
    N <- length(d)
    if(is.null(perm.data)){
      perm.data <- matrix(sample.int(N,size = N*permutations,replace = TRUE),nrow=N,ncol = permutations)
    }
    perm.data[]=d[perm.data]
    if(na.resample){isna=is.na(perm.data);if(sum(isna)>0){perm.data[isna]=sample(d[!is.na(d)],sum(isna))}}
    return(list(pool.sd=mean(sd.fun(perm.data)),pool.mean=mean(mean.fun(perm.data)),pool.se=sd(mean.fun(perm.data))))
  }
  pool.mix=if(ncol(data)==1){fun.mix(as.matrix(data[,1]),perm_mat,permutations)}else{sapply(1:ncol(data),function(i, ...) fun.mix(data[,i],perm_mat,permutations))}
  pool.mix <- df2array(as.matrix(pool.mix))
  colnames(pool.mix)=colnames(data)
  return(pool.mix)
}
"df2array" <-  function(df){
  return(matrix(unlist(df),nrow(df),ncol(df),dimnames = dimnames(df)))
}

##################
# functions
############################
"pool" <- function(data.mat,permutations=1000,
                   stat=c("mean","sd","se","ci"),sample.set=1,
                   na.action=c("resample","rm"),confidance=0.05,
                   report.fullset=FALSE, repool.SE_CI=FALSE,  replace=TRUE,MARGIN=2){
  pool(rnorm(1000))
  if(na.action[1]=="resample"|is.null(na.action)){na.resample=TRUE;na.rm=FALSE}else{na.resample=FALSE;na.rm=TRUE}
  data = as.matrix(data.mat)
  if(MARGIN == 1){data = t(data)}

  if(sum(stat%in%c("mean","se","sd","ci"))>0){
    pool.mixed <- pool.mean_sd_se(data, MARGIN=2, permutations=permutations, na.resample=na.resample, na.rm=na.rm)
    pool.ci = pool.CI.ofMEAN.Fun(data, MARGIN=2, permutations=1000, confidance=0.05, na.resample=na.resample, na.rm=na.rm)
    CI.expected=rbind(lower.expect=qt(confidance/2,permutations,0)*pool.mixed["pool.se",],upper.expect=qt(1-confidance/2,permutations,0)*pool.mixed["pool.se",])+pool.mixed["pool.mean",]
  }
  return( rbind(pool.mixed,pool.ci,CI.expected))
}
