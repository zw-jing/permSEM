find.best.fit.model<-function(fit,data=NULL,edge.pvalue.filter=0.1, MI.cut = 1,op=NULL,missing="listwise",change.ori.model.path=TRUE){

  extracted.model=extract.sem.model(fit)
  Model=initial.model=extracted.model$model.nolable
  model.table=extracted.model$model.table

  fit = lavaan(model = Model,data = data,missing=missing)
  I=fitMeasures(fit,c("chisq","df","cfi", "rmsea","rmsea.pvalue"))
  Stepwise.regression=matrix(0,1000,12);colnames(Stepwise.regression)=c("discard","chisq","df","pvalue","max.p","Model","added","chisq1","df1","pvalue1","max.p1","Model1")
  Stepwise.regression=as.data.frame(Stepwise.regression)
  Stepwise.regression[1,]=c(NA,I["chisq"], I["df"],1-pchisq(q = I["chisq"],df = I["df"],log.p = F),NA,NA,NA,NA,NA,NA,NA,NA)
  n=1
  adding="noadding";re.add=FALSE
  re2=NULL
  for(cycle.k in 1:200){
    n=n+1
    Est <- parameterEstimates(fit, ci = FALSE, standardized = TRUE)
    if(!change.ori.model.path){
      judge_retain=Est$pvalue<edge.pvalue.filter&!is.na(Est$pvalue)|(paste0(Est$lhs,Est$op,Est$rhs)%in%model.table$experession.nolable)
    }else{
      judge_retain=Est$pvalue<edge.pvalue.filter&!is.na(Est$pvalue)
    }

    sig.edge=Est[judge_retain,]
    sig.edge$label[sig.edge$label!=""]=paste(sig.edge$label[sig.edge$label!=""],sep="")
    re=paste(sig.edge$lhs,sig.edge$op,sig.edge$rhs,sep="")
    Model=paste0(re,collapse ="; ")

    if(sum(!judge_retain)){
      insig.edge=Est[!judge_retain,]
      insig.edge$label[insig.edge$label!=""]=paste(insig.edge$label[insig.edge$label!=""],sep="")
      re1=paste(insig.edge$lhs,insig.edge$op,insig.edge$rhs,sep="")
      if(sum(re1%in%re2)){
        re1=sort(unique(re1[re1%in%re2]))
        discard=paste0(re1,collapse ="; ")
      }else{
        discard=""
      }

    }else{
      discard=""
    }
    #print(paste("discarding : ",discard))

    fit = lavaan(model = Model,data = data,missing=missing)
    I=fitMeasures(fit,c("chisq","df","pvalue ", "cfi", "rmsea","rmsea.pvalue"))
    Stepwise.regression$discard[n]=discard
    Stepwise.regression$Model[n]=Model
    Stepwise.regression[n,c("chisq", "df", "pvalue", "max.p")]=c(I["chisq"], I["df"],1-pchisq(q = I["chisq"],df = I["df"],log.p = F),
                                                                 max(as.numeric(Stepwise.regression[n-1,"max.p"]),1-pchisq(q = I["chisq"],df = I["df"],log.p = F),na.rm=TRUE))

    print( Stepwise.regression[n,][,c("discard","chisq","df","pvalue","max.p")])

    MI <- modificationIndices(fit,sort.=TRUE,minimum.value= (MI.cut+0.0001),op = op)

    rownames(MI)=MI$experession.nolable=paste(MI$lhs,MI$op,MI$rhs,sep="")
    temp=model.table[model.table$op=="~",]
    temp2=model.table[model.table$op=="~~",]

    MI=MI[!MI$experession.nolable%in%
            c(paste0(temp$rhs,"~~",temp$lhs),
              paste0(temp$lhs,"~~",temp$rhs),
              paste0(temp$rhs,temp$op,temp$lhs),
              paste0(temp2$lhs,"~",temp2$rhs),
              paste0(temp2$rhs,"~",temp2$lhs)
            ),]


    #sbuset.MI1$experession=model.table[rownames(sbuset.MI1),]$experession
    if(nrow(MI)>0){
      for(j in 1:length(MI$experession.nolable)){
        adding=MI$experession.nolable[j]
        re.add=grepl(pattern = adding,x = discard )&Stepwise.regression$pvalue[n]<=Stepwise.regression$pvalue[n-1]
        re2=sort(unique(c(re,adding)))
        Model1=paste0(re2,collapse ="; ")
        if(!(Model1%in%Stepwise.regression$Model1)){break}
      }

    }else{
      break
    }
    #print(paste("adding : " , adding))
    Stepwise.regression$Model1[n]=Model=Model1

    fit = lavaan(model = Model,data = data,missing=missing)
    I=fitMeasures(fit,c("chisq","df","pvalue ", "cfi", "rmsea","rmsea.pvalue"))
    Stepwise.regression$added[n]=adding

    Stepwise.regression[n,c("chisq1", "df1", "pvalue1", "max.p1")]=c(I["chisq"], I["df"],1-pchisq(q = I["chisq"],df = I["df"],log.p = F),
                                                                     max(as.numeric(Stepwise.regression[n-1,"max.p1"]),1-pchisq(q = I["chisq"],df = I["df"],log.p = F),na.rm=TRUE))
    if(Model1%in%Stepwise.regression$Model1[1:(n-1)]){
      break
    }
    print( Stepwise.regression[n,][,c("added","chisq1","df1","pvalue1","max.p1")])
  }

  Stepwise.regression=Stepwise.regression[1:n,]
  #if(max(Stepwise.regression$max.p,na.rm = TRUE)>max(Stepwise.regression$max.p1,na.rm = TRUE)){
  best_Model=Stepwise.regression[2:n,]$Model[which(round(Stepwise.regression[2:n,]$pvalue,10)==round(max(Stepwise.regression[2:n,]$max.p,na.rm = TRUE),10))[1]]
  #}else{
  #best_Model=Stepwise.regression[2:n,]$Model1[which(round(Stepwise.regression[2:n,]$pvalue1,10)==round(max(Stepwise.regression[2:n,]$max.p1,na.rm = TRUE),10))[1]]
  #}
  bestfit=lavaan(best_Model, data = data,missing=missing)
  bestfit.model.table=extract.sem.model(bestfit)
  addpaths=bestfit.model.table$model.table[!(bestfit.model.table$model.table$experession.nolable%in%model.table$experession.nolable),]
  return(list(bestfit=bestfit,best.fit.Model=best_Model,addpaths=addpaths,Stepwise.regression=Stepwise.regression,initial.model=list(initial.model=initial.model,model.table=model.table),bestfit.model.table=bestfit.model.table$model.table))
}




