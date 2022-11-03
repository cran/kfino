## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(kfino)
library(dplyr)
library(foreach)
library(parallel)
library(doParallel)

## -----------------------------------------------------------------------------
data(lambs)
myIDE<-unique(lambs$IDE)

print(myIDE)

## ----error=TRUE---------------------------------------------------------------
param=list(m0=NULL,
             mm=NULL,
             pp=NULL,
             aa=0.001,
             expertMin=10,
             expertMax=45,
             sigma2_m0=1,
             sigma2_mm=0.05,
             sigma2_pp=5,
             K=15,
             seqp=seq(0.4,0.7,0.1))

t0 <- Sys.time()
resu1<-list()

for (i in seq_along(myIDE)){
  print(myIDE[i])
  tp.test<-filter(lambs,IDE == myIDE[i])
  print(dim(tp.test))
  resu1[[i]]<-kfino_fit(datain=tp.test,
              Tvar="dateNum",Yvar="Poids",
              param=param,
              doOptim=TRUE)
}
Sys.time() - t0

print(length(resu1))


## ----error=TRUE---------------------------------------------------------------

param=list(m0=NULL,
             mm=NULL,
             pp=NULL,
             aa=0.001,
             expertMin=10,
             expertMax=45,
             sigma2_m0=1,
             sigma2_mm=0.05,
             sigma2_pp=5,
             K=15,
             seqp=seq(0.4,0.7,0.1))

t0<-Sys.time()

simpleCall<-function(datain,Index,Tvar,Yvar,param){
  datain<-as.data.frame(datain)
  ici<-unique(datain[,"IDE"])
  tp.data<-datain[ datain[,"IDE"] == ici[Index],]

  tp.resu<-kfino::kfino_fit(datain=tp.data,
              Tvar=Tvar,Yvar=Yvar,
              param=param,
              doOptim=TRUE)
  return(tp.resu)
}

ncores<-parallel::detectCores()
myCluster<-parallel::makeCluster(ncores - 1)
doParallel::registerDoParallel(myCluster)

resu2<-foreach(i=seq_along(myIDE), .packages="kfino") %dopar% 
            simpleCall(datain=lambs,
                       Index=i,
                       Tvar="dateNum",
                       Yvar="Poids",
                       param=param)

parallel::stopCluster(myCluster)
Sys.time() - t0

print(length(resu2))


## -----------------------------------------------------------------------------
identical(resu1,resu2)

## -----------------------------------------------------------------------------
sessionInfo()

