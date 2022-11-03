test_that("kfino - ggplot produced correctly - quali", {
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

  p1<-kfino_plot(resuin=resu1,typeG="quali",
                 Tvar="Day",Yvar="Poids",Ident="IDE")
  # Output type
  expect_s3_class(p1$data,"data.frame")
  expect_s3_class(p1$layers[[1]]$geom,"GeomPoint")
  expect_s3_class(p1$layers[[2]]$geom,"GeomLine")
  expect_s3_class(p1$layers[[3]]$geom,"GeomLine")
  expect_s3_class(p1$layers[[4]]$geom,"GeomLine")

  p2<-kfino_plot(resuin=resu1,typeG="quanti",
                 Tvar="Day",Yvar="Poids",Ident="IDE")
  # Output type
  expect_s3_class(p2$data,"data.frame")
  expect_s3_class(p2$layers[[1]]$geom,"GeomPoint")
  expect_s3_class(p2$layers[[2]]$geom,"GeomLine")
  expect_s3_class(p2$layers[[3]]$geom,"GeomLine")
  expect_s3_class(p2$layers[[4]]$geom,"GeomLine")

  p3<-kfino_plot(resuin=resu1,typeG="prediction",
                 Tvar="Day",Yvar="Poids",Ident="IDE")
  # Output type
  expect_s3_class(p3$data,"data.frame")
  expect_s3_class(p3$layers[[1]]$geom,"GeomPoint")
  expect_s3_class(p3$layers[[2]]$geom,"GeomLine")
  expect_s3_class(p3$layers[[3]]$geom,"GeomLine")
  expect_s3_class(p3$layers[[4]]$geom,"GeomLine")

})
