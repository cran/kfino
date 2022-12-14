% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kfino.R
\name{kfino_fit}
\alias{kfino_fit}
\title{kfino_fit a function to detect outlier with a Kalman Filtering approach}
\usage{
kfino_fit(
  datain,
  Tvar,
  Yvar,
  param = NULL,
  doOptim = TRUE,
  method = "ML",
  threshold = 0.5,
  kappa = 10,
  kappaOpt = 7,
  verbose = FALSE
)
}
\arguments{
\item{datain}{an input data.frame of one time course to study (unique IDE)}

\item{Tvar}{char, time column name in the data.frame datain, a numeric vector
Tvar should be expressed as a proportion of day in seconds}

\item{Yvar}{char, name of the variable to predict in the data.frame datain}

\item{param}{list, a list of initialization parameters}

\item{doOptim}{logical, if TRUE optimization of the initial parameters,
default TRUE}

\item{method}{character, the method used to optimize the initial parameters:
Expectation-Maximization algorithm `"EM"` (faster) or Maximization
Likelihood `"ML"` (more robust), default `"ML"`}

\item{threshold}{numeric, threshold to qualify an observation as outlier
according to the label_pred, default 0.5}

\item{kappa}{numeric, truncation setting for likelihood optimization over 
initial parameters, default 10}

\item{kappaOpt}{numeric, truncation setting for the filtering and outlier 
detection step with optimized parameters, default 7}

\item{verbose}{write details if TRUE (optional), default FALSE.}
}
\value{
a S3 list with two data frames and a list of vectors of
kfino results

detectOutlier: The whole input data set with the detected outliers 
                     flagged and the prediction of the analyzed variable. 
                     the following columns are joined to the columns 
                     present in the input data set:
 \describe{
  \item{prediction}{the parameter of interest - Yvar - predicted}
  \item{label_pred}{the probability of the value being well predicted}
  \item{lwr}{lower bound of the confidence interval of the predicted value}
  \item{upr}{upper bound of the confidence interval of the predicted value}
  \item{flag}{flag of the value (OK value, KO value (outlier), OOR value
              (out of range values defined by the user in `kfino_fit` with 
              `expertMin`, `expertMax` input parameters). If 
              flag == OOR the 4 previous columns are set to NA.}
 }

PredictionOK: A subset of `detectOutlier` data set with the predictions 
        of the analyzed variable on possible values (OK and KO values)

kfino.results: kfino results (a list of vectors containing the 
        prediction of the analyzed variable, the probability to be an 
        outlier, the likelihood, the confidence interval of 
        the prediction and the flag of the data) on input parameters that 
        were optimized if the user chose this option
}
\description{
kfino_fit a function to detect outlier with a Kalman Filtering approach
}
\details{
The initialization parameter list `param` contains:
\describe{
 \item{mm}{(optional) numeric, target weight, NULL if the user wants to 
           optimize it}
 \item{pp}{(optional) numeric, probability to be correctly weighed, NULL if 
           the user wants to optimize it}
 \item{m0}{(optional) numeric, initial weight, NULL if the user wants to 
           optimize it}
 \item{aa}{numeric, rate of weight change, default 0.001 }
 \item{expertMin}{numeric, the minimal weight expected by the user}
 \item{expertMax}{numeric, the maximal weight expected by the user}
 \item{sigma2_m0}{numeric, variance of m0, default 1}
 \item{sigma2_mm}{numeric, variance of mm, related to the unit of Tvar,
       default 0.05}
 \item{sigma2_pp}{numeric, variance of pp, related to the unit of Yvar,
       default 5}
 \item{K}{numeric, a constant value in the outlier function (trapezium),
          by default K=5}
 \item{seqp}{numeric vector, sequence of pp probability to be correctly 
             weighted. default seq(0.5,0.7,0.1)}
}
It should be given by the user based on their knowledge of the animal or the 
data set. All parameters are compulsory except m0, mm and pp that can be
optimized by the algorithm. In the optimization step, those three parameters
are initialized according to the input data (between the expert
range) using quantile of the Y distribution (varying between 0.2 and 0.8 for
m0 and 0.5 for mm). pp is a sequence varying between 0.5 and 0.7. A
sub-sampling is performed to speed the algorithm if the number of possible
observations studied is greater than 500. Optimization is performed using
`"EM"` or `"ML"` method.
}
\examples{
data(spring1)
library(dplyr)

# --- With Optimization on initial parameters - ML method
t0 <- Sys.time()
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

resu1<-kfino_fit(datain=spring1,
              Tvar="dateNum",Yvar="Poids",
              doOptim=TRUE,method="ML",param=param1,
              verbose=TRUE)
Sys.time() - t0

# --- Without Optimization on initial parameters
t0 <- Sys.time()
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
              verbose=FALSE)
Sys.time() - t0
}
