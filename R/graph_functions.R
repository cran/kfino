#' kfino_plot a graphical function for the result of a kfino run
#'
#' @param resuin a list resulting of the kfino algorithm
#' @param typeG char, type of graphic, either detection of outliers (with
#'              qualitative or quantitative display) or prediction. must be
#'              "quanti" or "quali" or "prediction"
#' @param Tvar char, time variable in the data.frame datain
#' @param Yvar char, variable which was analysed in the data.frame datain
#' @param Ident char, column name of the individual id to be analyzed
#' @param title char, a graph title
#' @param labelX char, a label for x-axis
#' @param labelY char, a label for y-axis
#'
#' @details The produced graphic can be, according to typeG:
#' \describe{
#'  \item{quali}{This plot shows the detection of outliers with a qualitative 
#'       rule: OK values (black), KO values (outliers, purple) and OOR values 
#'       (out of range values defined by the user in `kfino_fit`, red) }
#'  \item{quanti}{This plot shows the detection of outliers with a quantitative 
#'       display using the calculated probability of the kfino algorithm}
#'  \item{prediction}{This plot shows the prediction of the analyzed variable 
#'        plus the OK values. Prediction corresponds to E[X_{t} | Y_{1...t}] 
#'        for each time point t. Between 2 time points, we used a simple 
#'        linear interpolation.}
#' }
#'
#' @importFrom ggplot2 aes aes_string ggplot geom_point geom_line
#' @importFrom ggplot2 scale_color_manual ggtitle xlab ylab
#' @return a ggplot2 graphic
#' @export
#'
#' @examples
#' data(spring1)
#' library(dplyr)
#'
#' print(colnames(spring1))
#'
#' # --- Without Optimisation on initial parameters
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
#'               doOptim=FALSE)
#'
#' # flags are qualitative
#' kfino_plot(resuin=resu2,typeG="quali",
#'             Tvar="Day",Yvar="Poids",Ident="IDE",
#'             title="kfino spring1",
#'             labelX="Time (day)",labelY="Weight (kg)")
#'
#' # flags are quantitative
#' kfino_plot(resuin=resu2,typeG="quanti",
#'             Tvar="Day",Yvar="Poids",Ident="IDE")
#'
#' # predictions on OK values
#' kfino_plot(resuin=resu2,typeG="prediction",
#'             Tvar="Day",Yvar="Poids",Ident="IDE")
kfino_plot<-function(resuin,
                     typeG,
                     Tvar,
                     Yvar,
                     Ident,
                     title=NULL,
                     labelX=NULL,
                     labelY=NULL){

  # Existence check
  if (is.null(resuin[[1]])) {
    stop("NULL object - No graph to provide. Please check your input object.")
  }
  if (!typeG %in% c("quanti","quali","prediction")) {
    stop("This type is not allowed.")
  }
  if (typeG %in% c("prediction","quanti") & is.null(resuin[[2]])) {
    stop("NULL object - No graph to provide. Please check your input object.")
  }

  # Some formatting
  tp<-as.data.frame(resuin[[1]])
  myIDE<-unique(tp[,Ident])

  if (is.null(title)){
    tp.title1<-paste0("kfino outlier detection - ",myIDE)
    tp.title2<-paste0("kfino prediction - ",myIDE)
  } else {
    tp.title1<-paste0(title," - ",myIDE)
    tp.title2<-paste0(title," - ",myIDE)
  }
  
  if (is.null(labelX)) labelX<-Tvar
  if (is.null(labelY)) labelY<-Yvar

  # graphics
  if (typeG == "quali"){
    if (!is.null(resuin[[2]])){
      g1<-ggplot(tp,aes_string(x=Tvar))+
          geom_point( aes_string(y=Yvar,color="flag")) +
          geom_line(data=tp[!is.na(tp$prediction),], aes(y=.data$prediction)) +
          geom_line(data=tp[!is.na(tp$prediction),], aes(y=.data$lwr),
                    color="green") +
          geom_line(data=tp[!is.na(tp$prediction),], aes(y=.data$upr),
                    color="green") +
          scale_color_manual(values =
                            c("KO"="purple", "OK" = "black", "OOR"="red")) +
          ggtitle(tp.title1) + xlab(labelX) + ylab(labelY)
    } else {
      g1<-ggplot(tp,aes_string(x=Tvar))+
          geom_point( aes_string(y=Yvar,color="flag")) +
          scale_color_manual(values =
                            c("KO"="purple", "OK" = "black", "OOR"="red")) +
          ggtitle(tp.title1) + xlab(labelX) + ylab(labelY)
    }
    return(g1)
  } else if (typeG == "quanti"){
    if (!is.null(resuin[[2]])){
      g1<-ggplot(tp,aes_string(x=Tvar))+
          geom_point( aes_string(y=Yvar,color="label_pred")) +
          geom_line(data=tp[!is.na(tp$prediction),], aes(y=.data$prediction)) +
          geom_line(data=tp[!is.na(tp$prediction),], aes(y=.data$lwr),
                  color="green") +
          geom_line(data=tp[!is.na(tp$prediction),], aes(y=.data$upr),
                  color="green") +
          ggtitle(tp.title1) + xlab(labelX) + ylab(labelY)
      return(g1)
    }
  } else if (typeG == "prediction"){
    tp2<-filter(tp,.data$flag == "OK")

    g1<-ggplot(tp2,aes_string(x=Tvar))+
        geom_point( aes_string(y=Yvar)) +
        geom_line(data=tp2, aes(y=.data$prediction)) +
        geom_line(data=tp2, aes(y=.data$lwr),color="green") +
        geom_line(data=tp2, aes(y=.data$upr),color="green") +
        ggtitle(tp.title2) + xlab(labelX) + ylab(labelY)
    return(g1)
  }
}

#-------------------- End of file -----------------------------------
