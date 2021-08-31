as.symbol.pvalue<-function(pvalue,
                           sig.level=c(0.001,0.01,0.05,0.1),
                           sig.symbol=c("***","**","*","."),
                           in.sig.symbol="ns",in.sig.out=TRUE,digits = 2,
                           align=c("none","left","right")){
  x=as.vector(as.numeric(pvalue[!is.na(pvalue)]))
  symbol.p=vector(length = length(pvalue));
  if(sum(is.na(pvalue))){symbol.p[is.na(pvalue)]="NA"}
  sym.p=vector(length = sum(!is.na(pvalue)))
  sig.symbol=sortby(sig.symbol,sig.level)
  sig.level=sort(sig.level)
  
  for(i in 1:length(sig.level)){
    logic.1=x<=sig.level[i]&x>=0
    if(sum(logic.1)){
      sym.p[logic.1]=sig.symbol[i]
      x[logic.1]=-1
    }
  }
  if(in.sig.out){
    sym.p[x>-1]=round(x[x>-1],digits)
  }else{
    sym.p[x>-1]=in.sig.symbol[1]
  }
  symbol.p[!is.na(pvalue)]=sym.p
  symbol.p=char.align(symbol.p, align =align[1], digits=digits )
  return(symbol.p)
}


