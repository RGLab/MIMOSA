context("MIMOSA fitting")

source("helper-init.R")


result<-MIMOSA(NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
     data=E, method='EM',
     subset=RefTreat%in%'Treatment'&ANTIGEN%in%'ENV',
     ref=ANTIGEN%in%'ENV'&RefTreat%in%'Reference')

expect_that(result,is_a("MIMOSAResultList"))
expect_that(names(result),equals("ENV"))
expect_that(length(result),equals(1))

context("getZ")

expect_that(getZ(result),is_a("matrix"))
expect_that(nrow(getZ(result)),equals(150))
expect_that(ncol(getZ(result)),equals(2))

context("getW")

W<-getW(result)
expect_that(W,is_a("data.frame"))
expect_that(dim(W),equals(c(2,1)))
expect_that(colnames(W),equals("ENV"))
expect_that(rownames(W),equals(c("w.nonresp","w.resp")))

context("countsTable")

P<-countsTable(result,proportion=TRUE)
C<-countsTable(result,proportion=FALSE)

expect_that(P,is_a("matrix"))
expect_that(ncol(P),equals(4))
expect_that(nrow(P),equals(150))

expect_that(C,is_a("matrix"))
expect_that(ncol(C),equals(4))
expect_that(nrow(C),equals(150))

expect_that(sum(C%%1),equals(0))
expect_that(sum(P%%1),equals(249))

context("volcanoPlot")

expect_that(volcanoPlot(result),throws_error())
expect_that(volcanoPlot(result,CYTNUM-CYTNUM_REF),is_a("ggplot"))

context("pData")
expect_that(pData(result),is_a("data.frame"))
expect_that(pData(result),is_a("data.table"))
expect_that(colnames(pData(result)),equals(c("UID","TCELLSUBSET","CYTOKINE","RefTreat","ANTIGEN")))
expect_that(nrow(pData(result)),equals(150))

context("rstan model")
library(rstan)
set.seed(100)
N=250
ps = c(rbeta(N/2,0.025*500,500),rbeta(N/2,0.01*500,500))
pu = rbeta(N,0.01*500,500)
ns = rbinom(N,50000,ps)
nu = rbinom(N,50000,pu)
Ns = rep(50000,50)
Nu = rep(50000,50)
# b = mean((ns/(Ns+ns))[ns/(Ns+ns)>nu/(Nu+nu)])
# a = mean(nu/(Nu+nu))

# d = MIMOSA:::simulate2(w2 = 0.5,A0=1,B0=10000,As=10,Bs=10000,NS=50000,N0=50000,obs=100)
d=data.frame(Ns,Nu,ns,nu)
b = with(d,mean((ns/(Ns))[ns/(Ns)>nu/(Nu)]))
a = with(d,mean(nu/(Nu)))


data=list(Ns=d$Ns,Nu=d$Nu,ns=d$ns,nu=d$nu,S=nrow(d),a=a/100,b=b/100)
# fit = stan(file="exec/mimosa.stan",model_name = "MIMOSA",data = data,cores = 4,control = list(adapt_delta=0.9),iter = 2000)
# plot(y=colMeans(extract(fit,pars="z")[[1]]),x=colMeans(extract(fit,pars="ps")[[1]] - extract(fit,pars=c("pu"))[[1]]))
# plot(y=colMeans(extract(fit,pars="z")[[1]]),x=with(d,ns/Ns-nu/Nu))
# plot(y=colMeans(extract(fit,pars="z")[[1]]),x=colMeans(extract(fit,pars="ps")[[1]]))
# plot(colMeans(extract(fit,pars="ps")[[1]])-colMeans(extract(fit,pars="pu")[[1]]),ns/Ns-nu/Nu)
# plot(colMeans(extract(fit,pars="ps")[[1]]),ns/Ns)
# plot(fit,pars=c("mu"))
# plot(fit,pars=c("phi_u","phi_s"))
# plot(fit,pars=c("proportion"))
modl = stan_model(system.file("exec/mimosa.stan",package="MIMOSA"))
set.seed(101)
testthat::expect_success(expr = {fit_opt = optimizing(modl,data=data,init=list(mu=c(0.01,0.025),phi_u=list(100),phi_s=list(100)),draws=100,as_vector=FALSE,algorithm="LBFGS",hessian=TRUE)})
table(fit_opt$par$z)
# plot(y=fit_opt$par$z,fit_opt$par$ps-fit_opt$par$pu)
# fit_opt$par$proportion
# fit_opt$par$mu
# fit_opt$par$phi_s
# fit_opt$par$phi_u
