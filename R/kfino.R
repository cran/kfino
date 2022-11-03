#' kfino_fit a function to detect outlier with a Kalman Filtering approach
#' @param datain an input data.frame of one time course to study (unique IDE)
#' @param Tvar char, time column name in the data.frame datain, a numeric vector
#'             Tvar should be expressed as a proportion of day in seconds
#' @param Yvar char, name of the variable to predict in the data.frame datain
#' @param param list, a list of initialization parameters
#' @param doOptim logical, if TRUE optimization of the initial parameters,
#'                default TRUE
#' @param method character, the method used to optimize the initial parameters:
#'               Expectation-Maximization algorithm `"EM"` (faster) or Maximization
#'               Likelihood `"ML"` (more robust), default `"ML"`
#' @param threshold numeric, threshold to qualify an observation as outlier
#'        according to the label_pred, default 0.5
#' @param kappa numeric, truncation setting for likelihood optimization over 
#'        initial parameters, default 10
#' @param kappaOpt numeric, truncation setting for the filtering and outlier 
#'        detection step with optimized parameters, default 7
#' @param verbose write details if TRUE (optional), default FALSE.
#'
#' @details The initialization parameter list `param` contains:
#' \describe{
#'  \item{mm}{(optional) numeric, target weight, NULL if the user wants to 
#'            optimize it}
#'  \item{pp}{(optional) numeric, probability to be correctly weighed, NULL if 
#'            the user wants to optimize it}
#'  \item{m0}{(optional) numeric, initial weight, NULL if the user wants to 
#'            optimize it}
#'  \item{aa}{numeric, rate of weight change, default 0.001 }
#'  \item{expertMin}{numeric, the minimal weight expected by the user}
#'  \item{expertMax}{numeric, the maximal weight expected by the user}
#'  \item{sigma2_m0}{numeric, variance of m0, default 1}
#'  \item{sigma2_mm}{numeric, variance of mm, related to the unit of Tvar,
#'        default 0.05}
#'  \item{sigma2_pp}{numeric, variance of pp, related to the unit of Yvar,
#'        default 5}
#'  \item{K}{numeric, a constant value in the outlier function (trapezium),
#'           by default K=5}
#'  \item{seqp}{numeric vector, sequence of pp probability to be correctly 
#'              weighted. default seq(0.5,0.7,0.1)}
#' }
#' It should be given by the user based on their knowledge of the animal or the 
#' data set. All parameters are compulsory except m0, mm and pp that can be
#' optimized by the algorithm. In the optimization step, those three parameters
#' are initialized according to the input data (between the expert
#' range) using quantile of the Y distribution (varying between 0.2 and 0.8 for
#' m0 and 0.5 for mm). pp is a sequence varying between 0.5 and 0.7. A
#' sub-sampling is performed to speed the algorithm if the number of possible
#' observations studied is greater than 500. Optimization is performed using
#' `"EM"` or `"ML"` method.
#'
#' @importFrom stats dnorm quantile na.omit
#' @importFrom dplyr mutate filter left_join arrange %>%
#' @importFrom dplyr .data if_else row_number select
#'
#' @return a S3 list with two data frames and a list of vectors of
#' kfino results
#' 
#' @return detectOutlier: The whole input data set with the detected outliers 
#'                      flagged and the prediction of the analyzed variable. 
#'                      the following columns are joined to the columns 
#'                      present in the input data set:
#'  \describe{
#'   \item{prediction}{the parameter of interest - Yvar - predicted}
#'   \item{label_pred}{the probability of the value being well predicted}
#'   \item{lwr}{lower bound of the confidence interval of the predicted value}
#'   \item{upr}{upper bound of the confidence interval of the predicted value}
#'   \item{flag}{flag of the value (OK value, KO value (outlier), OOR value
#'               (out of range values defined by the user in `kfino_fit` with 
#'               `expertMin`, `expertMax` input parameters). If 
#'               flag == OOR the 4 previous columns are set to NA.}
#'  }
#' @return PredictionOK: A subset of `detectOutlier` data set with the predictions 
#'         of the analyzed variable on possible values (OK and KO values)
#' @return kfino.results: kfino results (a list of vectors containing the 
#'         prediction of the analyzed variable, the probability to be an 
#'         outlier, the likelihood, the confidence interval of 
#'         the prediction and the flag of the data) on input parameters that 
#'         were optimized if the user chose this option
#'
#' @export
#' @examples
#' data(spring1)
#' library(dplyr)
#'
#' # --- With Optimization on initial parameters - ML method
#' t0 <- Sys.time()
#' param1<-list(m0=NULL,
#'              mm=NULL,
#'              pp=NULL,
#'              aa=0.001,
#'              expertMin=30,
#'              expertMax=75,
#'              sigma2_m0=1,
#'              sigma2_mm=0.05,
#'              sigma2_pp=5,
#'              K=2,
#'              seqp=seq(0.5,0.7,0.1))
#'
#' resu1<-kfino_fit(datain=spring1,
#'               Tvar="dateNum",Yvar="Poids",
#'               doOptim=TRUE,method="ML",param=param1,
#'               verbose=TRUE)
#' Sys.time() - t0
#'
#' # --- Without Optimization on initial parameters
#' t0 <- Sys.time()
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
#' resu2<-kfino_fit(datain=spring1,
#'               Tvar="dateNum",Yvar="Poids",
#'               param=param2,
#'               doOptim=FALSE,
#'               verbose=FALSE)
#' Sys.time() - t0
kfino_fit<-function(datain,Tvar,Yvar,
                    param=NULL,
                    doOptim=TRUE,method="ML",
                    threshold=0.5,kappa=10,kappaOpt=7,
                    verbose=FALSE){

  if( any(is.null(param[["expertMin"]]) |
          is.null(param[["expertMax"]])) )
    stop('You have to define expertMin and expertMax.')
  if( any(!is.numeric(param[["expertMin"]]) |
          !is.numeric(param[["expertMax"]])) )
    stop('expertMin and expertMax must be numeric.')
  if( !is.numeric(t(datain[,Yvar])))
    stop('Input parameter Yvar must contain numeric values.')
  if( !is.numeric(t(datain[,Tvar])))
    stop('Input parameter Tvar must contain numeric values.')

  #-------------------------------------------------------------
  # load objects
  m0<-param[["m0"]]
  mm<-param[["mm"]]
  pp<-param[["pp"]]
  if(is.null(param[["aa"]])){ aa<-0.001 } else { aa<-param[["aa"]] }
  expertMin<-param[["expertMin"]]
  expertMax<-param[["expertMax"]]
  if(is.null(param[["sigma2_m0"]])){
    sigma2_m0<-1
  } else {sigma2_m0<-param[["sigma2_m0"]]}
  if(is.null(param[["sigma2_mm"]])){
    sigma2_mm<-0.05
  } else {sigma2_mm<-param[["sigma2_mm"]]}
  if(is.null(param[["sigma2_pp"]])){
    sigma2_pp<-5
  } else {sigma2_pp<-param[["sigma2_pp"]]}
  if(is.null(param[["K"]])){ K<-5 } else { K<-param[["K"]] }
  if(is.null(param[["seqp"]])){
    seqp<-seq(0.5,0.7,1)
  } else {seqp<-param[["seqp"]]}
  #-------------------------------------------------------------

  # Flag impossible weight values (expert knowledge)
  # the algorithm works on the dataset with possible values
  # Create an artificial column numbering the rows for joining purpose
  # OK: 'good' value - OOR: out of range value - KO: outlier value
  datain<-as.data.frame(datain)
  datain<-datain %>%
    mutate(rowNum=row_number(),
           flag1=if_else(.data[[Yvar]] > expertMin &
                           .data[[Yvar]] <= expertMax,"OK","OOR"))
  tp.dt<-datain %>% filter(.data$flag1 == "OK")

  Y<-tp.dt[,Yvar]
  N<-nrow(tp.dt)
  Tps<-tp.dt[,Tvar]

  #WARNING WARNING: AU lieu de calculer L qui est arrondi à 0,
  # je calcule 10^N fois Lu et 10^N q. Au total pr chaque donnée
  # je multiplie par 100 mais comme c'est l'ordre de grandeur de doutlier()
  # ca ne me parait pas disproportionné.
  # Au lieu de 10 et 10 je fais simplement sqrt(expertMax - expertMin)
  scalingC=sqrt(expertMax - expertMin)

  #------------------------------------------------------------------------
  # Optimisation on Initial parameters or not
  #------------------------------------------------------------------------
  if (doOptim == TRUE){

    if (N > 500){
      # optim with sub-sampling
      if (verbose){
        print("-------:")
        print("Optimization of initial parameters ")
        print("with sub-sampling and ML method - result:")
      }
      bornem0=quantile(Y[1:N/4], probs = c(.2, .8))
      m0opt=quantile(Y[1:N/4], probs = c(.5))
      mmopt=quantile(Y[(3*N/4):N], probs = c(.5))

      if (verbose){
        cat("range m0: ",bornem0,"\n")
        cat("Initial m0opt: ",m0opt,"\n")
        cat("Initial mmopt: ",mmopt,"\n")
      }
      popt=0.5
      #--- Saving datain before sub-sampling
      YY=Y
      TpsTps=Tps
      NN=N

      Subechant=sort(sample(1:NN,50))
      Y=YY[Subechant]
      Tps=TpsTps[Subechant]
      N=50
      Vopt=utils_likelihood(list(m0=m0opt,
                                 mm=mmopt,
                                 pp=popt,
                                 aa=aa,
                                 expertMin=expertMin,
                                 expertMax=expertMax,
                                 sigma2_mm=sigma2_mm,
                                 sigma2_m0=sigma2_m0,
                                 sigma2_pp=sigma2_pp,
                                 K=K),
                            Y=Y,Tps=Tps,N=N,
                            scalingC=scalingC,kappaOpt=kappaOpt)

      for (m0 in seq(bornem0[1],bornem0[2],2) ){
        for (mm in seq((m0-5),(m0+20),2) ){
          for (p in seqp){
            # A voir si 50 sous-echantillons au hasard suffisent. Comme dans
            # Robbins Monroe, permet aussi de reduire l'impact de la troncature
            Subechant=sort(sample(1:NN,50))
            Y=YY[Subechant]
            Tps=TpsTps[Subechant]

            V=utils_likelihood(list(m0=m0,
                                    mm=mm,
                                    pp=p,
                                    aa=aa,
                                    expertMin=expertMin,
                                    expertMax=expertMax,
                                    sigma2_mm=sigma2_mm,
                                    sigma2_m0=sigma2_m0,
                                    sigma2_pp=sigma2_pp,
                                    K=K),
                                Y=Y,Tps=Tps,N=N,
                               scalingC=scalingC,kappaOpt=kappaOpt)
            if (V > Vopt){
              Vopt=V
              m0opt=m0
              mmopt=mm
              popt=p
            }
          }
        }
      }

      # Le fait de prendre des sous-echantillons remarquables
      # (ie 50 au lieu de 200) ne divise que par 2 le temps de calculs
      # (au lieu de 4 comme imaginé cela vient p.e de la fonction sample())
      Y=YY
      Tps=TpsTps
      N=NN
      if (verbose){
        print("Optimized parameters with ML method: ")
        cat("Optimized m0: ",m0opt,"\n")
        cat("Optimized mm: ",mmopt,"\n")
        cat("Optimized pp: ",popt,"\n")
        print("-------:")
      }

      resultat=utils_fit(param=list(mm=mmopt,
                                    pp=popt,
                                    m0=m0opt,
                                    aa=aa,
                                    expertMin=expertMin,
                                    expertMax=expertMax,
                                    sigma2_m0=sigma2_m0,
                                    sigma2_mm=sigma2_mm,
                                    sigma2_pp=sigma2_pp,
                                    K=K),
                         threshold=threshold,Y=Y,Tps=Tps,N=N,kappa=kappa)

    } else if (N > 50){
      # optimization without sub-sampling, 2 methods, EM or ML
      if (method == "EM"){
        if (verbose){
          print("-------:")
          print("Optimization of initial parameters with EM method - result:")
          print("no sub-sampling performed:")
        }
        bornem0=quantile(Y[1:N/2], probs = c(.2, .8))
        if (verbose)  cat("range m0: ",bornem0,"\n")

          #--- par dichotomie
          # borne basse
          N_etape_EM=10
          m0_tmp=m0_low=bornem0[1]
          m_tmp=m0_tmp
          p_tmp=0.5
          diff_m0=diff_mm=diff_p=20
          k=1
          while (diff_m0 > 0.5 && diff_p > 0.0001 && diff_mm > 2){
            Res_EM=utils_EM(param=list(mm=m_tmp,
                                     pp=p_tmp,
                                     m0=m0_tmp,
                                     aa=aa,
                                     expertMin=expertMin,
                                     expertMax=expertMax,
                                     sigma2_m0=sigma2_m0,
                                     sigma2_mm=sigma2_mm,
                                     sigma2_pp=sigma2_pp,
                                     K=K),
                          kappaOpt=kappaOpt, Y=Y,Tps=Tps,N=N,scalingC=scalingC)
            diff_m0=abs(m0_tmp - Res_EM$m0[[1]])
            diff_p=abs(p_tmp - Res_EM$pp)
            diff_mm=abs(m_tmp - Res_EM$mm[[1]])
            k<-k+1
            
            m0_tmp=Res_EM$m0[[1]]
            m_tmp=Res_EM$mm[[1]]
            p_tmp=Res_EM$pp
            if (k==N_etape_EM) break
          }
          Vopt_low=Res_EM$likelihood
          m0opt_low<-Res_EM$m0[[1]]
          mmopt_low<-Res_EM$mm[[1]]
          popt_low<-Res_EM$pp
          
          # borne haute
          N_etape_EM=10
          m0_tmp=m0_up=bornem0[2]
          m_tmp=m0_tmp
          p_tmp=0.5
          diff_m0=diff_mm=diff_p=20
          k=1
          while (diff_m0 > 0.5 && diff_p > 0.0001 && diff_mm > 2){
            Res_EM=utils_EM(param=list(mm=m_tmp,
                                     pp=p_tmp,
                                     m0=m0_tmp,
                                     aa=aa,
                                     expertMin=expertMin,
                                     expertMax=expertMax,
                                     sigma2_m0=sigma2_m0,
                                     sigma2_mm=sigma2_mm,
                                     sigma2_pp=sigma2_pp,
                                     K=K),
                          kappaOpt=kappaOpt, Y=Y,Tps=Tps,N=N,scalingC=scalingC)
            diff_m0=abs(m0_tmp - Res_EM$m0[[1]])
            diff_p=abs(p_tmp - Res_EM$pp)
            diff_mm=abs(m_tmp - Res_EM$mm[[1]])
            k<-k+1
            
            m0_tmp=Res_EM$m0[[1]]
            m_tmp=Res_EM$mm[[1]]
            p_tmp=Res_EM$pp
            if (k==N_etape_EM) break
          }
          Vopt_up=Res_EM$likelihood
          m0opt_up<-Res_EM$m0[[1]]
          mmopt_up<-Res_EM$mm[[1]]
          popt_up<-Res_EM$pp
          
          diff_m0range<-abs(m0opt_up - m0opt_low)
          diff_pprange<-abs(popt_up - popt_low)
          diff_mmrange<-abs(mmopt_up - mmopt_low)
          
          # test quasi equality
          while(diff_m0range > 0.5 && diff_pprange > 0.0001 && diff_mmrange > 2){
            m0_tmp=m0_med=(m0_up+m0_low)/2
            m_tmp=m0_tmp
            p_tmp=0.5
            diff_m0=diff_mm=diff_p=20
            k=1
            while (diff_m0 > 0.5 && diff_p > 0.0001 && diff_mm > 2){
               print(k)
              Res_EM=utils_EM(param=list(mm=m_tmp,
                                       pp=p_tmp,
                                       m0=m0_tmp,
                                       aa=aa,
                                       expertMin=expertMin,
                                       expertMax=expertMax,
                                       sigma2_m0=sigma2_m0,
                                       sigma2_mm=sigma2_mm,
                                       sigma2_pp=sigma2_pp,
                                       K=K),
                                    kappaOpt=kappaOpt, Y=Y,Tps=Tps,N=N,
                                    scalingC=scalingC)
              diff_m0=abs(m0_tmp - Res_EM$m0[[1]])
              diff_p=abs(p_tmp - Res_EM$pp)
              diff_mm=abs(m_tmp - Res_EM$mm[[1]])
              k<-k+1
              
              m0_tmp=Res_EM$m0[[1]]
              m_tmp=Res_EM$mm[[1]]
              p_tmp=Res_EM$pp
              if (k==N_etape_EM) break
            }
            
            if (Vopt_up < Vopt_low){
              Vopt_up<-Res_EM$likelihood
              m0opt_up<-Res_EM$m0[[1]]
              mmopt_up<-Res_EM$mm[[1]]
              popt_up<-Res_EM$pp
              m0_up<-m0_med
            } else {
              Vopt_low<-Res_EM$likelihood
              m0opt_low<-Res_EM$m0[[1]]
              mmopt_low<-Res_EM$mm[[1]]
              popt_low<-Res_EM$pp
              m0_low<-m0_med
            }
            
            diff_m0range<-abs(m0opt_up - m0opt_low)
            diff_pprange<-abs(popt_up - popt_low)
            
          }
          
          # challenger
          m0opt=quantile(Y[1:N/4], probs = c(.5))
          mmopt=quantile(Y[(3*N/4):N], probs = c(.5))
          popt=0.5
          
          Vopt=utils_EM(param=list(m0=m0opt,
                                 mm=mmopt,
                                 pp=popt,
                                 aa=aa,
                                 expertMin=expertMin,
                                 expertMax=expertMax,
                                 sigma2_mm=sigma2_mm,
                                 sigma2_m0=sigma2_m0,
                                 sigma2_pp=sigma2_pp,
                                 K=K),
                      Y=Y,Tps=Tps,N=N,scalingC=scalingC,
                      kappaOpt=kappaOpt)$likelihood
          
          if (Vopt_low > Vopt){
            m0opt<-m0opt_low
            mmopt<-mmopt_low
            popt<-popt_low
          }
          
          if (verbose){
            print("Optimized parameters with EM method: ")
            cat("Optimized m0: ",m0opt,"\n")
            cat("Optimized mm: ",mmopt,"\n")
            cat("Optimized pp: ",popt,"\n")
            print("-------:")
          }
          resultat=utils_fit(param=list(mm=mmopt,
                                        pp=popt,
                                        m0=m0opt,
                                        aa=aa,
                                        expertMin=expertMin,
                                        expertMax=expertMax,
                                        sigma2_m0=sigma2_m0,
                                        sigma2_mm=sigma2_mm,
                                        sigma2_pp=sigma2_pp,
                                        K=K),
                             threshold=threshold,Y=Y,Tps=Tps,N=N,kappa=kappa)
          
      } else if (method == "ML"){
        if (verbose){
          print("-------:")
          print("Optimization of initial parameters with ML method - result:")
          print("no sub-sampling performed:")
        }
        bornem0=quantile(Y[1:N/4], probs = c(.2, .8))
        m0opt=quantile(Y[1:N/4], probs = c(.5))
        mmopt=quantile(Y[(3*N/4):N], probs = c(.5))

        if (verbose){
          cat("range m0: ",bornem0,"\n")
          cat("initial m0opt: ",m0opt,"\n")
          cat("initial mmopt: ",mmopt,"\n")
        }
        popt=0.5

        Vopt=utils_likelihood(list(m0=m0opt,
                                  mm=mmopt,
                                  pp=popt,
                                  aa=aa,
                                  expertMin=expertMin,
                                  expertMax=expertMax,
                                  sigma2_mm=sigma2_mm,
                                  sigma2_m0=sigma2_m0,
                                  sigma2_pp=sigma2_pp,
                                  K=K),
                              Y=Y,Tps=Tps,N=N,
                              scalingC=scalingC,kappaOpt=kappaOpt)
        for (m0 in seq(bornem0[1],bornem0[2],2) ){
          for (mm in seq((m0-5),(m0+20),2) ){
            for (p in seqp){
              V=utils_likelihood(list(m0=m0,
                                      mm=mm,
                                      pp=p,
                                      aa=aa,
                                      expertMin=expertMin,
                                      expertMax=expertMax,
                                      sigma2_mm=sigma2_mm,
                                      sigma2_m0=sigma2_m0,
                                      sigma2_pp=sigma2_pp,
                                      K=K),
                                 Y=Y,Tps=Tps,N=N,
                                 scalingC=scalingC,kappaOpt=kappaOpt)

              if (V > Vopt){
                Vopt=V
                m0opt=m0
                mmopt=mm
                popt=p
              }
            }
          }
        }

        if (verbose){
          print("Optimized parameters: ")
          cat("Optimized m0: ",m0opt,"\n")
          cat("Optimized mm: ",mmopt,"\n")
          cat("Optimized pp: ",popt,"\n")
          print("-------:")
        }
        resultat=utils_fit(param=list(mm=mmopt,
                                      pp=popt,
                                      m0=m0opt,
                                      aa=aa,
                                      expertMin=expertMin,
                                      expertMax=expertMax,
                                      sigma2_m0=sigma2_m0,
                                      sigma2_mm=sigma2_mm,
                                      sigma2_pp=sigma2_pp,
                                      K=K),
                           threshold=threshold,Y=Y,Tps=Tps,N=N,kappa=kappa)
      }

    } else {
      if (N > 5) {
        # Not enough data - no optim
        if (is.null(param[["m0"]])){
          X<-c(trunc(quantile(Y[1:N/4], probs = .2)), 0.5,
               round(quantile(Y[(3*N/4):N], probs = 0.8)))
        } else {
          X<-c(m0,pp,mm)
        }
        
        if (verbose){
          print("-------:")
          print("Optimization of initial parameters - result:")
          print("Not enough data => No optimization performed:")
          print("Used parameters: ")
          print(X)
          print("-------:")
        }
        resultat=utils_fit(param=list(m0=X[[1]],
                                      pp=X[[2]],
                                      mm=X[[3]],
                                      aa=aa,
                                      expertMin=expertMin,
                                      expertMax=expertMax,
                                      sigma2_m0=sigma2_m0,
                                      sigma2_mm=sigma2_mm,
                                      sigma2_pp=sigma2_pp,
                                      K=K),
                           threshold=threshold,Y=Y,Tps=Tps,N=N,kappa=kappa)
      } else {
        warning("Not enough data between expert knowledge. The algorithm is not performed.")
        resultat<-NULL
      }
    }
  } else {
    # Pas d'optimisation et test si N petit - si trop petit on ne fait rien
    if (N > 5) {
      # No optimisation on initial parameters
      if (is.null(param[["m0"]])){
        X<-c(trunc(quantile(Y[1:N/4], probs = .2)), 0.5,
             round(quantile(Y[(3*N/4):N], probs = 0.8)))
      } else {
        X<-c(m0,pp,mm)
      }
      if (verbose){
        print("-------:")
        print("No optimization of initial parameters:")
        print("Used parameters: ")
        print(X)
      }
      resultat=utils_fit(param=list(m0=X[[1]],
                                    pp=X[[2]],
                                    mm=X[[3]],
                                    aa=aa,
                                    expertMin=expertMin,
                                    expertMax=expertMax,
                                    sigma2_m0=sigma2_m0,
                                    sigma2_mm=sigma2_mm,
                                    sigma2_pp=sigma2_pp,
                                    K=K),
                         threshold=threshold,Y=Y,Tps=Tps,N=N,kappa=kappa)
    } else {
      warning("Not enough data between expert knowledge. The algorithm is not performed.")
      resultat<-NULL
    }
  }

  #---------------------------------------------------------------------------
  # Formatting output results
  #---------------------------------------------------------------------------
  # If resultat NULL then create output with datain and 2 NULL objects
  # else create output with detectoutlier, PredictionOK and kfino.results
  # useful for the kfino_plot() function
  #---------------------------------------------------------------------------
  if (is.null(resultat)){
    dt.out<-datain %>% 
              mutate(flag=.data$flag1) %>%
              select(-.data$flag1)
    dt.pred<-NULL
    resultat<-NULL

    mylist<-list(dt.out,dt.pred,resultat)
    names(mylist)<-c("detectOutlier","PredictionOK","kfino.results")
    class(mylist) = c("kfino")
    return(invisible(mylist))
  } else {
    prediction=na.omit(resultat$prediction)
    label_pred=round(na.omit(resultat$label),2)
    lwr=na.omit(resultat$lwr)
    upr=na.omit(resultat$upr)
    flag=na.omit(resultat$flag)

    dt.pred=cbind.data.frame(select(tp.dt,.data$rowNum),
                             prediction,label_pred,lwr,upr,flag)

    # join dt.pred with the initial datain containing all the observations
    dt.out<-left_join(datain,dt.pred,by="rowNum")
    dt.out<-mutate(dt.out,flag=if_else(.data$flag1 == "OOR",
                                       .data$flag1, .data$flag))
    dt.out<-arrange(dt.out,.data$rowNum)
    dt.out<-select(dt.out,-.data$flag1)

    #--------------------------------------
    # return a S3 list object
    #--------------------------------------
    # 1. a whole dataset with the detected outliers flagged and prediction
    # 2. a dataset with the prediction on possible values
    # 3. optimization results (a list of vectors)
    mylist<-list(dt.out,dt.pred,resultat)
    names(mylist)<-c("detectOutlier","PredictionOK","kfino.results")
    class(mylist) = c("kfino")
    return(invisible(mylist))
  }
}



#-------------------------------------- End of file --------------------------
