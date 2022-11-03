## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(kfino)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
data(spring1)

# Dimension of this dataset
dim(spring1)

head(spring1)

## ----error=TRUE---------------------------------------------------------------
# --- Without Optimisation on parameters
param2<-list(m0=41,
             mm=45,
             pp=0.5,
             aa=0.001,
             expertMin=30,
             expertMax=75,
             sigma2_m0=1,
             sigma2_mm=0.05,
             sigma2_pp=5,
             K=2,
             seqp=seq(0.5,0.7,0.1))

resu2<-kfino_fit(datain=spring1,
              Tvar="dateNum",Yvar="Poids",
              param=param2,
              doOptim=FALSE,
              verbose=TRUE)     

## -----------------------------------------------------------------------------
# structure of detectOutlier data set
str(resu2$detectOutlier)

# head of PredictionOK data set
head(resu2$PredictionOK)

# structure of kfino.results list
str(resu2$kfino.results)

## -----------------------------------------------------------------------------
# flags are qualitative
kfino_plot(resuin=resu2,typeG="quali",
            Tvar="Day",Yvar="Poids",Ident="IDE")
            
# flags are quantitative
kfino_plot(resuin=resu2,typeG="quanti",
            Tvar="Day",Yvar="Poids",Ident="IDE")

## -----------------------------------------------------------------------------
# --- With Optimisation on parameters
param1<-list(m0=NULL,
             mm=NULL,
             pp=NULL,
             aa=0.001,
             expertMin=30,
             expertMax=75,
             sigma2_m0=1,
             sigma2_mm=0.05,
             sigma2_pp=5,
             K=2,
             seqp=seq(0.5,0.7,0.1))


## ----error=TRUE---------------------------------------------------------------

resu1<-kfino_fit(datain=spring1,
              Tvar="dateNum",Yvar="Poids",
              param=param1,
              doOptim=TRUE,
              method="ML",
              verbose=TRUE)  

# flags are qualitative
kfino_plot(resuin=resu1,typeG="quali",
            Tvar="Day",Yvar="Poids",Ident="IDE")
            
# flags are quantitative
kfino_plot(resuin=resu1,typeG="quanti",
            Tvar="Day",Yvar="Poids",Ident="IDE")

## ----error=TRUE---------------------------------------------------------------
kfino_plot(resuin=resu1,typeG="prediction",
            Tvar="Day",Yvar="Poids",Ident="IDE")


## ----error=TRUE---------------------------------------------------------------
resu1b<-kfino_fit(datain=spring1,
              Tvar="dateNum",Yvar="Poids",
              param=param1,
              doOptim=TRUE,
              method="EM",
              verbose=TRUE)  

# flags are qualitative
kfino_plot(resuin=resu1b,typeG="quali",
            Tvar="Day",Yvar="Poids",Ident="IDE")
            
# flags are quantitative
kfino_plot(resuin=resu1b,typeG="quanti",
            Tvar="Day",Yvar="Poids",Ident="IDE")

kfino_plot(resuin=resu1b,typeG="prediction",
            Tvar="Day",Yvar="Poids",Ident="IDE")


## -----------------------------------------------------------------------------
data(merinos1)

# Dimension of this dataset
dim(merinos1)

head(merinos1)

## ----error=TRUE---------------------------------------------------------------
# --- With Optimisation on parameters
param3<-list(m0=NULL,
             mm=NULL,
             pp=NULL,
             aa=0.001,
             expertMin=10,
             expertMax=45,
             sigma2_m0=1,
             sigma2_mm=0.05,
             sigma2_pp=5,
             K=2,
             seqp=seq(0.5,0.7,0.1))

resu3<-kfino_fit(datain=merinos1,
              Tvar="dateNum",Yvar="Poids",
              doOptim=TRUE,param=param3,
              verbose=TRUE)      

# flags are qualitative
kfino_plot(resuin=resu3,typeG="quali",
            Tvar="Day",Yvar="Poids",Ident="IDE")
            
# flags are quantitative
kfino_plot(resuin=resu3,typeG="quanti",
            Tvar="Day",Yvar="Poids",Ident="IDE")

## ----error=TRUE---------------------------------------------------------------
kfino_plot(resuin=resu3,typeG="prediction",
            Tvar="Day",Yvar="Poids",Ident="IDE")


## ----session,echo=FALSE,message=FALSE, warning=FALSE--------------------------
  sessionInfo()

