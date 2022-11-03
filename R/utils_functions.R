#-------------------------------------------------------------------
# utils_functions.R: some useful functions for kfino method
# doutlier()
# utils_fit()
# utils_likelihood()
# utils_EM()
#-------------------------------------------------------------------

#' doutlier defines an outlier distribution (Surface of a
#' trapezium) and uses input parameters given in the main function kfino_fit()
#'
#' @param y numeric, point
#' @param K numeric, constant value
#' @param expertMin numeric, the minimal weight expected by the user
#' @param expertMax numeric, the maximal weight expected by the user
#'
#' @details this function is used to calculate an outlier distribution
#'          following a trapezium shape. 
#'  \eqn{y \mapsto \text{doutlier}(y,K,\text{expertMin},\text{expertMax})}
#'          is the probability density function on 
#'  \eqn{[\text{expertMin},\text{expertMax}]} which is linear and verifies
#'  \eqn{\text{doutlier}(\text{expertMax},K,\text{expertMin},\text{expertMax}) 
#'  =K*\text{doutlier}(\text{expertMin},K,\text{expertMin},\text{expertMax}).}
#'  In particular, when $K=1$ this corresponds to the uniform distribution.
#'  
#' @return a numeric value
#' @export
#'
#' @examples
#' doutlier(2,5,10,45)
doutlier<-function(y,
                      K,
                      expertMin,
                      expertMax){
  2/((K+1)*(expertMax - expertMin)) +
    (2*(K-1)/(K+1))*((y - expertMin)/((expertMax - expertMin)^2))
}


#--------------------------------------------------------------------------
#' utils_fit a fonction running the kfino algorithm to filter data and 
#'   detect outliers under the knowledge of all parameters
#'
#' @param param list, see initial parameter list in \code{kfino_fit}
#' @param threshold numeric, threshold for confidence interval, default 0.5
#' @param kappa numeric, truncation setting for likelihood optimization, 
#'        default 10
#' @param Y character, name of the numeric variable to predict in the  
#'        data.frame datain
#' @param Tps character, time column name in the data.frame datain, a 
#'            numeric vector.
#'            Tvar can be expressed as a proportion of day in seconds
#' @param N numeric, length of the numeric vector of Y values
#'
#' @details utils_fit is a tool function used in the main \code{kfino_fit} 
#' function. It uses the same input parameter list than the main function.
#' @return a list
#' \describe{
#'  \item{prediction}{vector, the prediction of weights}
#'  \item{label}{vector, probability to be an outlier}
#'  \item{likelihood}{numeric, the calculated likelihood}
#'  \item{lwr}{vector of lower bound confidence interval of the prediction }
#'  \item{upr}{vector of upper bound confidence interval of the prediction }
#'  \item{flag}{char, is an outlier or not}
#' }
#' @export
#'
#' @examples
#' set.seed(1234)
#' Y<-rnorm(n=10,mean=50,4)
#' Tps<-seq(1,10)
#' N=10
#' param2<-list(m0=41,
#'              mm=45,
#'              pp=0.5,
#'              aa=0.001,
#'              expertMin=30,
#'              expertMax=75,
#'              sigma2_m0=1,
#'              sigma2_mm=0.05,
#'              sigma2_pp=5,
#'              K=2,
#'             seqp=seq(0.5,0.7,0.1))
#' print(Y)
#' utils_fit(param=param2,threshold=0.5,kappa=10,Y=Y,Tps=Tps,N=N)
utils_fit<-function(param,threshold,kappa=10,Y,Tps,N){
  # load objects
  mm=param[["mm"]]
  pp=param[["pp"]]
  m0=param[["m0"]]
  aa=param[["aa"]]
  expertMin=param[["expertMin"]]
  expertMax=param[["expertMax"]]
  sigma2_m0=param[["sigma2_m0"]]
  sigma2_mm=param[["sigma2_mm"]]
  sigma2_pp=param[["sigma2_pp"]]
  K=param[["K"]]
  # paramètre de troncature
  kappa <-kappa

  # initialisation (1.1.1)
  #--------------------
  m1= (sigma2_pp*m0 + Y[1]*sigma2_m0)/(sigma2_m0+sigma2_pp)
  sigma1=(sigma2_m0*sigma2_pp)/(sigma2_m0+sigma2_pp)

  l0<- doutlier(Y[1],K,expertMin,expertMax)
  loinorm1<-dnorm(Y[1],m0, sqrt(sigma2_m0+sigma2_pp))

  p0= ((1-pp)*l0) / (pp*loinorm1 + (1-pp)*l0)
  p1=1-p0

  m=c(m0,m1)
  p=c(p0,p1)
  sigma2=c(sigma2_m0,sigma1)
  L=c(l0,loinorm1)

  Xmean=c(p0*m0 + p1*m1)
  Pmean=c(p1)
  q=c(1,1)
  sigma2_new=c(sigma2_m0)

  #iteration (1.1.2)
  #-----------------------
  #// Pour l'instant, je fais comme si kappa<N-1 mettre un if sinon
  #before truncation
  for (k in 1:(kappa-1)){
    # je crée le vecteur des nouveaux m, ##pour plus tard: au lieu de les
    # faire vide, je prend direct m_{u0}
    mnew=rep(0,2^(k+1))
    sigma2new=rep(0 ,2^(k+1))
    pnew=rep(0 ,2^(k+1))
    Lnew=rep(0 ,2^(k+1))
    # CR: renormalization constant that intervenes in the denominator of the pu
    qnew=rep(0 ,2^(k+1))
    diffTps<-Tps[k+1] - Tps[k]
    #--- numérateur de pu0
    tpbeta<-doutlier(Y[k+1],K,expertMin,expertMax)

    pnew[1:(2^k)]=p[1:(2^k)]*(1-pp)*tpbeta
    Lnew[1:(2^k)]=L[1:(2^k)]*tpbeta
    mnew[1:(2^k)]= m[1:(2^k)]*exp(-aa*diffTps) + mm*(1-exp(-aa*diffTps)) #m_u0

    sigma2new[1:(2^k)] = sigma2[1:(2^k)]*exp(-2*aa*(diffTps)) +
      (1-exp(-2*aa*(diffTps)))*sigma2_mm/(2*aa)
    qnew[1:(2^k)] <- q[1:(2^k)]*(1-pp)
    qnew[(1+2^k):2^(k+1)] <- q[1:(2^k)]*pp
    sommevar=sigma2new[1:(2^k)] + sigma2_pp
    tpnorm<-dnorm(Y[k+1], mnew[1:(2^k)], sqrt(sommevar))
    pnew[(1+2^k):2^(k+1)]=p[1:(2^k)]*pp*tpnorm
    Lnew[(1+2^k):2^(k+1)]=L[1:(2^k)]*tpnorm
    mnew[(1+2^k):2^(k+1)] = (sigma2new[1:(2^k)]*Y[k+1] +
                               mnew[1:(2^k)]*sigma2_pp)/sommevar
    sigma2new[(1+2^k):2^(k+1)] = sigma2new[1:(2^k)] * sigma2_pp/sommevar

    CR<-sum(pnew)
    Proba1<-sum(pnew[(1+2^k):2^(k+1)])
    SommeMasse<-sum(pnew*mnew)

    m=mnew
    sigma2=sigma2new
    p=pnew/CR
    Xmean=c(Xmean,SommeMasse/CR)
    Pmean=c(Pmean,Proba1/CR)
    L=Lnew
    q=qnew
    sigma2_new<-c(sigma2_new,sum(p*sigma2) + sum(p*m^2) - (sum(p*m))^2)
  }

  # after truncation
  #----------------------
  for (k in kappa:(N-1)){
    # Initialization of the new m vector
    # TODO: instead of 0, perhaps set directly at m_{u0}
    mnew=rep(0,2^(kappa+1))
    sigma2new=rep(0 ,2^(kappa+1))
    pnew=rep(0 ,2^(kappa+1))
    Lnew=rep(0 ,2^(kappa+1))
    # CR: renormalization constant that intervenes in the denominator of the pu
    qnew=rep(0 ,2^(kappa+1))
    diffTps<-Tps[k+1] - Tps[k]

    #--- pu0 numerator
    tpbeta<-doutlier(Y[k+1],K,expertMin,expertMax)

    pnew[1:(2^kappa)]=p[1:(2^kappa)]*(1-pp)*tpbeta
    Lnew[1:(2^kappa)]=L[1:(2^kappa)]*tpbeta
    mnew[1:(2^kappa)]= m[1:(2^kappa)]*exp(-aa*diffTps) + mm*(1-exp(-aa*diffTps)) #m_u0
    sigma2new[1:(2^kappa)] = sigma2[1:(2^kappa)]*exp(-2*aa*(diffTps)) +
      (1-exp(-2*aa*(diffTps)))*sigma2_mm/(2*aa)
    qnew[1:(2^kappa)]=q[1:(2^kappa)]*(1-pp)
    qnew[(1+2^kappa):2^(kappa+1)]=q[1:(2^kappa)]*pp
    sommevar=sigma2new[1:(2^kappa)] + sigma2_pp
    tpnorm<-dnorm(Y[k+1], mnew[1:(2^kappa)], sqrt(sommevar))
    pnew[(1+2^kappa):2^(kappa+1)]=p[1:(2^kappa)]*pp*tpnorm
    Lnew[(1+2^kappa):2^(kappa+1)]=L[1:(2^kappa)]*tpnorm
    mnew[(1+2^kappa):2^(kappa+1)] = (sigma2new[1:(2^kappa)]*Y[k+1] +
                                       mnew[1:(2^kappa)]*sigma2_pp)/sommevar #m_u1
    sigma2new[(1+2^kappa):2^(kappa+1)] = sigma2new[1:(2^kappa)] * sigma2_pp/sommevar

    CR<-sum(pnew)
    Proba1<-sum(pnew[(1+2^kappa):2^(kappa+1)])
    SommeMasse<-sum(pnew*mnew)

    selection=order(pnew, decreasing=T)[1:2^kappa]

    m=mnew[selection]
    sigma2=sigma2new[selection]
    p=pnew[selection]/sum(pnew[selection])
    Xmean=c(Xmean,SommeMasse/CR)
    Pmean=c(Pmean,Proba1/CR)
    L=Lnew[selection]
    q=qnew[selection]
    sigma2_new<-c(sigma2_new,sum(p*sigma2) + sum(p*m^2) - (sum(p*m))^2)
  }

  #----------------------------------------------
  # Ajout Intervalle de confiance des estimations
  # Variance (see section 1.2.1)
  # IC est approximé...
  #----------------------------------------------
  lwr <- Xmean - 1.96*sqrt(sigma2_new)
  upr <- Xmean + 1.96*sqrt(sigma2_new)
  # flag
  flag<-if_else(Pmean > threshold,"OK","KO")

  #--- Calcul des derniers éléments
  Vraisemblance=L%*%q

  resultat=list("prediction"=Xmean,
                "label"=Pmean,
                "likelihood"=Vraisemblance,
                "lwr"=lwr,
                "upr"=upr,
                "flag"=flag)
  return(resultat)
}

#------------------------------------------------------------------------
#' utils_likelihood a function to calculate a likelihood on initial parameters 
#' optimized by a grid search
#' 
#' @param param list, see initial parameter list in \code{kfino_fit}
#' @param kappaOpt numeric, truncation setting for initial parameters' 
#'        optimization, default 7
#' @param Y character, name of the numeric variable to predict in the  
#'        data.frame datain
#' @param Tps character, time column name in the data.frame datain, a 
#'            numeric vector.
#'            Tvar can be expressed as a proportion of day in seconds
#' @param N numeric, length of the numeric vector of Y values
#' @param scalingC numeric, scaling constant. To be changed if the function is 
#'  not able to calculate the likelihood because the number of data is large
#'
#' @details utils_likelihood is a tool function used in the main 
#' \code{kfino_fit} function. It uses the same input parameter list than 
#' the main function.
#' @return a likelihood
#' @keywords internal
#' @export
#'
#' @examples
#' set.seed(1234)
#' Y<-rnorm(n=10,mean=50,4)
#' Tps<-seq(1,10)
#' N=10
#' param2<-list(m0=41,
#'              mm=45,
#'              pp=0.5,
#'              aa=0.001,
#'              expertMin=30,
#'              expertMax=75,
#'              sigma2_m0=1,
#'              sigma2_mm=0.05,
#'              sigma2_pp=5,
#'              K=2,
#'              seqp=seq(0.5,0.7,0.1))
#' print(Y)
#' utils_likelihood(param=param2,kappaOpt=7,Y=Y,Tps=Tps,N=N,scalingC=6)
utils_likelihood<-function(param,kappaOpt=7,Y,Tps,N,scalingC){
  # load objects
  mm=param[["mm"]]
  pp=param[["pp"]]
  m0=param[["m0"]]
  aa=param[["aa"]]
  expertMin=param[["expertMin"]]
  expertMax=param[["expertMax"]]
  sigma2_m0<-param[["sigma2_m0"]]
  sigma2_mm<-param[["sigma2_mm"]]
  sigma2_pp<-param[["sigma2_pp"]]
  K<-param[["K"]]

  #---- truncation parameter
  # Here, kappa is set to 7, the probabilities are < 0.01 instead of 0.001
  # with kappa=10 and so computing time is divided by 10.
  kappa=kappaOpt

  # initialisation (1.1.1)
  #--------------------
  m1= (sigma2_pp*m0 + Y[1]*sigma2_m0)/(sigma2_m0+sigma2_pp)
  sigma1=(sigma2_m0*sigma2_pp)/(sigma2_m0+sigma2_pp)

  l0<- doutlier(Y[1],K,expertMin,expertMax)
  loinorm1<-dnorm(Y[1],m0, sqrt(sigma2_m0+sigma2_pp))

  p0= ((1-pp)*l0) / (pp*loinorm1 + (1-pp)*l0)
  p1=1-p0

  m=c(m0,m1)
  p=c(p0,p1)
  sigma2=c(sigma2_m0,sigma1)
  L=scalingC*c(l0,loinorm1) # increase it with the scaling constant
  q=c(1,1) 

  #iteration (1.1.2)
  #-----------------------
  # before truncation
  for (k in 1:(kappa-1)){
    mnew=rep(0,2^(k+1))
    sigma2new=rep(0 ,2^(k+1))
    pnew=rep(0 ,2^(k+1))
    Lnew=rep(0 ,2^(k+1))
    # CR: renormalization constant that intervenes in the denominator of the pu
    qnew=rep(0 ,2^(k+1))
    diffTps<-Tps[k+1] - Tps[k]
    #--- numerator of pu0
    tpbeta<-doutlier(Y[k+1],K,expertMin,expertMax)
    pnew[1:(2^k)]=p[1:(2^k)]*(1-pp)*tpbeta
    Lnew[1:(2^k)]=L[1:(2^k)]*tpbeta
    mnew[1:(2^k)]= m[1:(2^k)]*exp(-aa*diffTps) + mm*(1-exp(-aa*diffTps)) #m_u0
    sigma2new[1:(2^k)] = sigma2[1:(2^k)]*exp(-2*aa*(diffTps)) +
      (1-exp(-2*aa*(diffTps)))*sigma2_mm/(2*aa)
    qnew[1:(2^k)] <- q[1:(2^k)]*(1-pp)
    qnew[(1+2^k):2^(k+1)] <- q[1:(2^k)]*pp
    sommevar=sigma2new[1:(2^k)] + sigma2_pp
    tpnorm<-dnorm(Y[k+1], mnew[1:(2^k)], sqrt(sommevar))
    pnew[(1+2^k):2^(k+1)]=p[1:(2^k)]*pp*tpnorm
    Lnew[(1+2^k):2^(k+1)]=L[1:(2^k)]*tpnorm
    mnew[(1+2^k):2^(k+1)] = (sigma2new[1:(2^k)]*Y[k+1] + mnew[1:(2^k)]*sigma2_pp)/sommevar
    sigma2new[(1+2^k):2^(k+1)] = sigma2new[1:(2^k)] * sigma2_pp/sommevar

    m=mnew
    sigma2=sigma2new
    p=pnew/sum(pnew)

    L=scalingC*Lnew 
    q=scalingC*qnew
  }

  #  after truncation
  #----------------------
  for (k in kappa:(N-1)){
    # creation of the 'new m' vector
    mnew=rep(0,2^(kappa+1))
    sigma2new=rep(0 ,2^(kappa+1))
    pnew=rep(0 ,2^(kappa+1))
    Lnew=rep(0 ,2^(kappa+1))
    # CR: renormalization constant that intervenes in the denominator of the pu
    qnew=rep(0 ,2^(kappa+1))
    diffTps<-Tps[k+1] - Tps[k]
    #--- numerator of pu0
    tpbeta<-doutlier(Y[k+1],K,expertMin,expertMax)
    pnew[1:(2^kappa)]=p[1:(2^kappa)]*(1-pp)*tpbeta
    Lnew[1:(2^kappa)]=L[1:(2^kappa)]*tpbeta
    mnew[1:(2^kappa)]= m[1:(2^kappa)]*exp(-aa*diffTps) + mm*(1-exp(-aa*diffTps)) #m_u0
    sigma2new[1:(2^kappa)] = sigma2[1:(2^kappa)]*exp(-2*aa*(diffTps)) +
      (1-exp(-2*aa*(diffTps)))*sigma2_mm/(2*aa)
    qnew[1:(2^kappa)]=q[1:(2^kappa)]*(1-pp)
    qnew[(1+2^kappa):2^(kappa+1)]=q[1:(2^kappa)]*pp
    sommevar=sigma2new[1:(2^kappa)] + sigma2_pp
    tpnorm<-dnorm(Y[k+1], mnew[1:(2^kappa)], sqrt(sommevar))
    pnew[(1+2^kappa):2^(kappa+1)]=p[1:(2^kappa)]*pp*tpnorm
    Lnew[(1+2^kappa):2^(kappa+1)]=L[1:(2^kappa)]*tpnorm
    mnew[(1+2^kappa):2^(kappa+1)] = (sigma2new[1:(2^kappa)]*Y[k+1] +
                                       mnew[1:(2^kappa)]*sigma2_pp)/sommevar #m_u1
    sigma2new[(1+2^kappa):2^(kappa+1)] = sigma2new[1:(2^kappa)] * sigma2_pp/sommevar

    selection=order(pnew, decreasing=T)[1:2^kappa]

    m=mnew[selection]
    sigma2=sigma2new[selection]
    p=pnew[selection]/sum(pnew[selection])
    L=scalingC*Lnew[selection] 
    q=scalingC*qnew[selection]
  }

  Vraisemblance=L%*%q

  return(Vraisemblance)
}

#------------------------------------------------------------------------
#' utils_EM a function to estimate the parameters `m_0` , `mm`, `pp` through 
#' an Expectation-Maximization (EM) method
#'
#' @param param list, see initial parameter list in \code{kfino_fit}
#' @param kappaOpt numeric, truncation setting for initial parameters' 
#'        optimization, default 7
#' @param Y character, name of the numeric variable to predict in the  
#'        data.frame datain
#' @param Tps character, time column name in the data.frame datain, a 
#'            numeric vector.
#'            Tvar can be expressed as a proportion of day in seconds
#' @param N numeric, length of the numeric vector of Y values
#' @param scalingC numeric, scaling constant. To be changed if the function is 
#'  not able to calculate the likelihood because the number of data is large
#'
#' @details utils_EM is a tool function used in the main \code{kfino_fit} 
#' function. It uses the same input parameter list than the main function.
#' @return a list:
#' \describe{
#'  \item{m0}{numeric, optimized m0}
#'  \item{mm}{numeric, optimized mm}
#'  \item{pp}{numeric, optimized pp}
#'  \item{likelihood}{numeric, the calculated likelihood}
#' }
#' @export
#'
#' @examples
#' set.seed(1234)
#' Y<-rnorm(n=10,mean=50,4)
#' Tps<-seq(1,10)
#' N=10
#' param2<-list(m0=41,
#'              mm=45,
#'              pp=0.5,
#'              aa=0.001,
#'              expertMin=30,
#'              expertMax=75,
#'              sigma2_m0=1,
#'              sigma2_mm=0.05,
#'              sigma2_pp=5,
#'              K=2,
#'              seqp=seq(0.5,0.7,0.1))
#' print(Y)
#' utils_EM(param=param2,kappaOpt=7,Y=Y,Tps=Tps,N=N,scalingC=6)
utils_EM<-function(param,kappaOpt,Y,Tps,N,scalingC){
  # load objects
  mm<-param[["mm"]]
  pp<-param[["pp"]]
  m0<-param[["m0"]]
  aa<-param[["aa"]]
  expertMin<-param[["expertMin"]]
  expertMax<-param[["expertMax"]]
  sigma2_m0<-param[["sigma2_m0"]]
  sigma2_mm<-param[["sigma2_mm"]]
  sigma2_pp<-param[["sigma2_pp"]]
  K<-param[["K"]]
  
  #---- truncation parameter
  # Here, kappa is set to 7, the probabilities are < 0.01 instead of 0.001
  # with kappa=10 and so computing time is divided by 10.
  kappa=kappaOpt
  
  # initialisation (1.1.1)
  #--------------------
  m1= (sigma2_pp*m0 + Y[1]*sigma2_m0)/(sigma2_m0+sigma2_pp)
  sigma1=(sigma2_m0*sigma2_pp)/(sigma2_m0+sigma2_pp)
  
  l0<- doutlier(Y[1],K,expertMin,expertMax)
  loinorm1<-dnorm(Y[1],m0, sqrt(sigma2_m0+sigma2_pp))
  
  p0= ((1-pp)*l0) / (pp*loinorm1 + (1-pp)*l0)
  p1=1-p0
  
  m=c(m0,m1)
  p=c(p0,p1)
  sigma2=c(sigma2_m0,sigma1)
  L=scalingC*c(l0,loinorm1) 
  q=c(1,1) 
  
  a=c(1,sigma2_pp/(sigma2_m0+sigma2_pp))
  b=c(0,0)
  c=c(0, Y[1]*sigma2_m0/(sigma2_m0+sigma2_pp))
  matrice_A=matrix(a*a/(sigma2 + sigma2_pp), ncol=1, nrow=2, byrow=TRUE)
  matrice_B=matrix(b*b/(sigma2 + sigma2_pp), ncol=1, nrow=2, byrow=TRUE)
  matrice_C=matrix(a*b/(sigma2 + sigma2_pp), ncol=1, nrow=2, byrow=TRUE)
  matrice_Ya=matrix((a*(Y[1]-c)/(sigma2 + sigma2_pp)), ncol=1, nrow=2, byrow=TRUE)
  matrice_Yb=matrix((b*(Y[1]-c)/(sigma2 + sigma2_pp)), ncol=1, nrow=2, byrow=TRUE)
  
  Z=matrix(c(0,1), ncol=1, nrow=2, byrow=TRUE)
  #iteration (1.1.2)
  #-----------------------
  # kappa < N-1 for now. add an if condition if necessary
  # before truncation
  for (k in 1:(kappa-1)){
    mnew=rep(0,2^(k+1))
    sigma2new=rep(0 ,2^(k+1))
    pnew=rep(0 ,2^(k+1))
    Lnew=rep(0 ,2^(k+1))
    anew=rep(0 ,2^(k+1))
    bnew=rep(0 ,2^(k+1))
    cnew=rep(0 ,2^(k+1))
    
    # CR: renormalization constant that intervenes in the denominator of the pu
    qnew=rep(0 ,2^(k+1))
    diffTps<-Tps[k+1] - Tps[k]
    #--- numérateur de pu0
    tpbeta<-doutlier(Y[k+1],K,expertMin,expertMax)
    pnew[1:(2^k)]=p[1:(2^k)]*(1-pp)*tpbeta
    Lnew[1:(2^k)]=L[1:(2^k)]*tpbeta
    mnew[1:(2^k)]= m[1:(2^k)]*exp(-aa*diffTps) + mm*(1-exp(-aa*diffTps)) #m_u0
    sigma2new[1:(2^k)] = sigma2[1:(2^k)]*exp(-2*aa*(diffTps)) +
      (1-exp(-2*aa*(diffTps)))*sigma2_mm/(2*aa)
    qnew[1:(2^k)] <- q[1:(2^k)]*(1-pp)
    qnew[(1+2^k):2^(k+1)] <- q[1:(2^k)]*pp
    # new parameters
    anew[1:(2^k)]= a[1:(2^k)]*exp(-aa*diffTps)
    bnew[1:(2^k)]= b[1:(2^k)]*exp(-aa*diffTps) + (1-exp(-aa*diffTps))
    cnew[1:(2^k)]= c[1:(2^k)]*exp(-aa*diffTps)
    
    # continuation
    sommevar=sigma2new[1:(2^k)] + sigma2_pp
    tpnorm<-dnorm(Y[k+1], mnew[1:(2^k)], sqrt(sommevar))
    pnew[(1+2^k):2^(k+1)]=p[1:(2^k)]*pp*tpnorm
    Lnew[(1+2^k):2^(k+1)]=L[1:(2^k)]*tpnorm
    mnew[(1+2^k):2^(k+1)] = (sigma2new[1:(2^k)]*Y[k+1] + mnew[1:(2^k)]*sigma2_pp)/sommevar
    sigma2new[(1+2^k):2^(k+1)] = sigma2new[1:(2^k)] * sigma2_pp/sommevar
    # new parameters
    anew[(1+2^k):2^(k+1)]= a[1:(2^k)]*sigma2_pp/sommevar
    bnew[(1+2^k):2^(k+1)]= b[1:(2^k)]*sigma2_pp/sommevar
    cnew[(1+2^k):2^(k+1)]= c[1:(2^k)]*sigma2_pp/sommevar + Y[k+1]*sigma2new[1:(2^k)]/sommevar
    
    m=mnew
    sigma2=sigma2new
    p=pnew/sum(pnew)
    # new parameters
    a=anew
    b=bnew
    c=cnew
    
    matrice_A=cbind(rbind(matrice_A,matrice_A), a*a/(sigma2 + sigma2_pp)) 
    matrice_B=cbind(rbind(matrice_B,matrice_B), b*b/(sigma2 + sigma2_pp)) 
    matrice_C=cbind(rbind(matrice_C,matrice_C), a*b/(sigma2 + sigma2_pp)) 
    matrice_Ya=cbind(rbind(matrice_Ya,matrice_Ya), a*(Y[k+1]-c)/(sigma2 + sigma2_pp)) 
    matrice_Yb=cbind(rbind(matrice_Yb,matrice_Yb), b*(Y[k+1]-c)/(sigma2 + sigma2_pp)) 
    
    Znew=matrix(rep(0,(k+1)*2^(k+1)), ncol=k+1, nrow=2^(k+1), byrow=TRUE)
    Znew[(1:2^k),]=cbind(Z, rep(0,2^k))
    Znew[(2^k+1):(2^(k+1)),]=cbind(Z, rep(1,2^k))
    Z=Znew
    
    L=scalingC*Lnew # fois 2 pr le grandir
    q=scalingC*qnew
  }
  
  # after truncation
  #----------------------
  for (k in kappa:(N-1)){
    # Initialisation of the new m vectors
    mnew=rep(0,2^(kappa+1))
    sigma2new=rep(0 ,2^(kappa+1))
    pnew=rep(0 ,2^(kappa+1))
    Lnew=rep(0 ,2^(kappa+1))
    anew=rep(0 ,2^(kappa+1))
    bnew=rep(0 ,2^(kappa+1))
    cnew=rep(0 ,2^(kappa+1))
    # CR: renormalization constant that intervenes in the denominator of the pu
    qnew=rep(0 ,2^(kappa+1))
    diffTps<-Tps[k+1] - Tps[k]
    #--- pu0 numerator
    tpbeta<-doutlier(Y[k+1],K,expertMin,expertMax)
    pnew[1:(2^kappa)]=p[1:(2^kappa)]*(1-pp)*tpbeta
    Lnew[1:(2^kappa)]=L[1:(2^kappa)]*tpbeta
    # m_uO
    mnew[1:(2^kappa)]= m[1:(2^kappa)]*exp(-aa*diffTps) + mm*(1-exp(-aa*diffTps))
    sigma2new[1:(2^kappa)] = sigma2[1:(2^kappa)]*exp(-2*aa*(diffTps)) +
      (1-exp(-2*aa*(diffTps)))*sigma2_mm/(2*aa)
    qnew[1:(2^kappa)]=q[1:(2^kappa)]*(1-pp)
    qnew[(1+2^kappa):2^(kappa+1)]=q[1:(2^kappa)]*pp
    # new parameters
    anew[1:(2^kappa)]= a[1:(2^kappa)]*exp(-aa*diffTps)
    bnew[1:(2^kappa)]= b[1:(2^kappa)]*exp(-aa*diffTps) + (1-exp(-aa*diffTps))
    cnew[1:(2^kappa)]= c[1:(2^kappa)]*exp(-aa*diffTps)
    sommevar=sigma2new[1:(2^kappa)] + sigma2_pp
    tpnorm<-dnorm(Y[k+1], mnew[1:(2^kappa)], sqrt(sommevar))
    pnew[(1+2^kappa):2^(kappa+1)]=p[1:(2^kappa)]*pp*tpnorm
    Lnew[(1+2^kappa):2^(kappa+1)]=L[1:(2^kappa)]*tpnorm
    # m_u1
    mnew[(1+2^kappa):2^(kappa+1)] = (sigma2new[1:(2^kappa)]*Y[k+1] +
                                     mnew[1:(2^kappa)]*sigma2_pp)/sommevar
    sigma2new[(1+2^kappa):2^(kappa+1)]=sigma2new[1:(2^kappa)] * sigma2_pp/sommevar
    
    # new parameters
    anew[(1+2^kappa):2^(kappa+1)]= a[1:(2^kappa)]*sigma2_pp/sommevar
    bnew[(1+2^kappa):2^(kappa+1)]= b[1:(2^kappa)]*sigma2_pp/sommevar
    cnew[(1+2^kappa):2^(kappa+1)]= c[1:(2^kappa)]*sigma2_pp/sommevar + 
                                   Y[k+1]*sigma2new[1:(2^kappa)]/sommevar
    
    selection=order(pnew, decreasing=T)[1:2^kappa]
    
    m=mnew[selection]
    sigma2=sigma2new[selection]
    p=pnew[selection]/sum(pnew[selection])
    L=scalingC*Lnew[selection] 
    q=scalingC*qnew[selection]
    
    Znew=matrix(rep(0,(k+1)*2^(kappa+1)), ncol=k+1, nrow=2^(kappa+1), byrow=TRUE)
    Znew[1:2^kappa,]=cbind(Z, rep(0,2^kappa))
    Znew[((1+2^kappa):2^(kappa+1)),]=cbind(Z, rep(1,2^kappa))
    
    Z=Znew[selection,]
    
    a=anew[selection]
    b=bnew[selection]
    c=cnew[selection]
    matrice_A=(cbind(rbind(matrice_A,matrice_A), 
                     anew*anew/(sigma2new + sigma2_pp)))[selection,]
    matrice_B=(cbind(rbind(matrice_B,matrice_B), 
                     bnew*bnew/(sigma2new + sigma2_pp)))[selection,]
    matrice_C=(cbind(rbind(matrice_C,matrice_C), 
                     anew*bnew/(sigma2new + sigma2_pp)))[selection,]
    matrice_Ya=(cbind(rbind(matrice_Ya,matrice_Ya), 
                      anew*(Y[k+1]-cnew)/(sigma2new + sigma2_pp)))[selection,]
    matrice_Yb=(cbind(rbind(matrice_Yb,matrice_Yb), 
                      bnew*(Y[k+1]-cnew)/(sigma2new + sigma2_pp)))[selection,]
  }
  
  Vraisemblance=L%*%q
  
  A=(Z[,1]%*%p)/(sigma2_m0+sigma2_pp) + sum((p%*%(matrice_A*Z))[1:(N-1)])
  Ya=(Z[,1]%*%p)*Y[1]/(sigma2_m0+sigma2_pp) + sum((p%*%(matrice_Ya*Z))[1:(N-1)])
  B= sum((p%*%(matrice_B*Z))[1:(N-1)])
  C= sum((p%*%(matrice_C*Z))[1:(N-1)])
  Yb= sum((p%*%(matrice_Yb*Z))[1:(N-1)])
  
  newm0=(C*Yb-B*Ya)/(C^2-A*B)
  newmm=(C*Ya-A*Yb)/(C^2-A*B)
  ppnew=max(sum(p%*%Z)/N , 0.5)
  
  # Outputs
  resultat=list("m0"=newm0,
                "mm"=newmm, 
                "pp"=ppnew, 
                "likelihood"=Vraisemblance)
  return(resultat)
}


#------------------------ End of file -----------------------------------
