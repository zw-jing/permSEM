\encoding{UTF-8}
\name{lavaan.path.stat}
\alias{lavaan.path.stat}
\title{testing  paths in SEMs}
\description{(1) Testing the significanes direct/indirect/total effects of factor (or grouped factors) on other factor (or grouped factors)in SEMs paths [sem.effects.stat(), lavaan.path.stat()]; (2) Helping the users to build modelings [add.variance(), delete.variables.in.model()] more efficiently; (3) Helping the users to find the better/best fitting modelings based on the initial modelings [find.best.fit.model()] (4) See more detailed descriptions with "?lavaan.path.stat"}
\usage{
lavaan.path.stat(fit,data, effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                 group.rhs=NULL, group.lhs=NULL, standardized = TRUE,
                 via=NULL,
                 bootstrap=FALSE,  permutation=FALSE, reps=1000, return.boot=TRUE,
                 missing="listwise", ind.in.dist=NULL, dist.list=NULL, boot.perm.simultaneous = FALSE,
                 strata=NULL, exclude.paths=NULL, exclude.pars=NULL)
}
\arguments{
  \item{fit}{a lavaan fit structure}
  \item{data}{dataframe with all the variables used in SEM}
  \item{effect.type}{indicating significances of one or multiple effect of c("total.effects", "direct.effects", and "indirect.effects") to be tested}
  \item{group.rhs}{rhs to be grouped together. e.g., New.grouped.rhs = "rh1+rh2+rh3"}
  \item{group.lhs}{lhs to be grouped together. e.g., New.grouped.lhs = "lh1+lh2+lh3"}
  \item{via}{indicating the factors, through which indirect effects of lhs on rhs were realized}
  \item{bootstrap}{whether perferform bootstrap to test significances of effects}
  \item{reps}{randomizions for bootstrap and/or permutation test}
  \item{return.boot}{wheter return bootstraped (or obserbed + permuted null) effects}
  \item{permutation}{whether perferform a mantel like permutation test based on distance matrices}
  \item{dist.list}{a list of distance or similarity matrices}
  \item{ind.in.dist}{index of subset data of distances used for SEM}
  \item{boot.perm.simultaneous}{whether perferform bootstrap (for calculating unbiased mean) and permutation (for testing significance) simultaneously or not}

  \item{strata}{within which group permutation (random sampling) are to be performed}
  \item{exclude.paths}{paths (e.g., "rh1~lh1+lh2") not to be included for estimitions of indirect effects}
  \item{exclude.pars}{factors (e.g., c("f1","f2")") not to be included for estimitions of indirect effects}
}

\author{
Zhongwang Jing
}

\examples{
data(C_N_coupling)
library(lavaan)
library(semPlot)
library(vegan)
library(permSEM)
df=C_N_coupling[C_N_coupling$crop=="Rice" & C_N_coupling$N.fert=="FN",]
df[,c(8:ncol(df))]=vegan::decostand(df[,c(8:ncol(df))],MARGIN = 2,method = "standardize",na.rm = T)
names(df)

df$TNs=df$TNs-df$DNs
df$TNf=df$TNf-df$DNf

###################  buid initial sem model
#基于已知经验或理论构建初始模型,已知的路径都要包含在模型里（参考：Oberski, D. lavaan.survey: An R Package for Complex Survey Analysis of Structural Equation Models. Journal of Statistical Software 57, 1-27 (2014).）
model.1 <-
  '
TNf+TNs+DNf+DNs~LUC+Nitro+Rhizo
SOC12+SOC13+DOC~LUC+Nitro+Rhizo
TNf+TNs+DNf+DNs~SOC12+SOC13+DOC
DNs~TNs
DNf~TNf
narG+nosZ+ AOB+AOA+nirK+ MBC+Bac~LUC+Nitro+SOC13+DOC+Rhizo
TNf+TNs+DNf+DNs~narG+nosZ+ AOB+AOA+nirK+ MBC+Bac
'
#LFLF.PD
model.1=add.variance(model = model.1,covariance = c("Bac","nirK","narG","AOA","AOB","nosZ"))
saturated_fit=lavaan(model = model.1 ,data = as.data.frame(vegan::decostand(C_N_coupling[,5:ncol(C_N_coupling)],MARGIN = 2,method = "standardize",na.rm = T)))

########################
model.1=delete.variables.in.model(saturated_fit,delete.pars = c("LUC", "Nitro" ,
                                                                "SOC12","nosZ",
                                                                "AOA","AOB","Rhizo",
                                                                "DN","TN","MBC"))

fit=lavaan(model = model.1, data = df, ordered = NULL,
           sampling.weights   = NULL,
           sample.cov = NULL, sample.mean = NULL, sample.th = NULL,
           sample.nobs = NULL,
           group = NULL, cluster = NULL)
summary(fit)
#fitMeasures(fit)
fitMeasures(fit,c("chisq","df","pvalue ", "cfi", "rmsea","rmsea.pvalue"))
#modificationIndices(fit,sort.=TRUE,minimum.value= 1)
semPaths(fit,layout = "circle2",what = "est",whatLabels = "std.all",style="ram",residuals=TRUE,nCharNodes=6,
         sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)

#semPaths(fit,layout = "circle2",whatLabels = "std.all")

#查看直接间接效应（自己写的代码，供参考）
effects.stat=sem.effects.stat(fit = fit,data = df)
#effects.stat$effects.summary
effects.stat1=lavaan.path.stat(fit = fit,data = df)

effects.stat$std.effects.summary$total.effects
effects.stat1$std.effects.summary$total.effects

effects.stat$std.effects.summary$direct.effects
effects.stat1$std.effects.summary$direct.effects

effects.stat$std.effects.summary$indirect.effects
effects.stat1$std.effects.summary$indirect.effects

#在初始模型的基础上做逐步回归，找最优模型 （自己写的代码，普适性可能会较差，可能会出现bug；最初目的是为模型的构建提供参考）
MD=find.best.fit.model(fit,data = df,op = "~~",change.ori.model.path = F,edge.pvalue.filter = 0.1)
#change.ori.model.path : 是(TRUE)否(FALSE)删去初始模型中不显著的路径
#edge.pvalue.filter: 路径的p值阈值，大于该阈值考虑将其删去
#op : 添加单向路径（A~B: B影响A）；添加相关（covariance）（A~~B: A与B的残差相关）；NULL: 自动选择添加影响（~）或相关（~~）路径

bestfit = lavaan(model = MD$best.fit.Model, data = df)
summary(bestfit)
semPaths(bestfit,layout = "circle2",what = "est",whatLabels = "std.all",style="ram",residuals=TRUE,nCharNodes=6,
         sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)
print(sem.effects.stat(fit = bestfit,data = df)$std.effects.summary)

eff = lavaan.path.stat(bestfit,data=df,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                       group.rhs=c(#Rhizo="Rhizo",
                         #LUC="LUC",
                         Carbon="SOC13+DOC",
                         Bio.process="MBC+narG+nirK+nosZ+AOB+AOA+Bac",
                         TNs="TNs",TNf="TNf"),
                       group.lhs=c(TNs="TNs",TNf="TNf",
                                   DNs="DNs",DNf="DNf",DN="DN",
                                   Bio.process="MBC+narG+nirK+nosZ+AOB+AOA+Bac"),
                       via="MBC+narG+nirK+nosZ+AOB+AOA+Bac+TNf+TNs",
                       rh=NULL,lh=NULL,standardized = TRUE,
                       bootstrap=TRUE,reps=10,return.boot=TRUE,boot.depth = NULL,permutation = F,
                       ind.in.dist = NULL,dist.list = NULL,boot.perm.simultaneous = F)

eff$std.effects.summary


}
