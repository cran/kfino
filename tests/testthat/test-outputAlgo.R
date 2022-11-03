test_that("kfino - no optimisation - computed correctly", {
  data(spring1)
  param1<-list(m0=41,
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
  resu1<-kfino_fit(datain=spring1,
                    Tvar="dateNum",Yvar="Poids",
                    param=param1,
                    doOptim=FALSE)

  # Output type
  expect_equal(length(resu1),3)
  expect_s3_class(resu1[[1]],"data.frame")
  expect_s3_class(resu1[[2]],"data.frame")
  expect_type(resu1[[3]],"list")
  expect_equal(length(resu1[[3]]),6)

  # expected result
  expect_equal(as.vector(resu1[[3]]$likelihood),1.245621e-150)
})


test_that("kfino - optimisation - computed correctly", {
  data(spring1)
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
                    param=param1,
                    doOptim=TRUE)

  # Output type
  expect_equal(length(resu1),3)
  expect_s3_class(resu1[[1]],"data.frame")
  expect_s3_class(resu1[[2]],"data.frame")
  expect_type(resu1[[3]],"list")
  expect_equal(length(resu1[[3]]),6)

  # expected result
  expect_equal(as.vector(resu1[[3]]$likelihood),1.035855e-135)
})

test_that("utils_likelihood - likelihood - computed correctly", {
  set.seed(1234)
  Y<-rnorm(n=10,mean=50,4)
  Tps<-seq(1,10)
  N=10
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
  resu<-utils_likelihood(param=param2,Y=Y,Tps=Tps,N=N,scalingC=6,kappaOpt=7)

  expect_equal(as.vector(resu),0.00490257592)

})

