context("MIMOSA fitting")

source("helper-init.R")


result<- suppressWarnings(
  MIMOSA(
    NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
    data=E, method='EM',
    subset=RefTreat%in%'Treatment'&ANTIGEN%in%'ENV',
    ref=ANTIGEN%in%'ENV'&RefTreat%in%'Reference'
    )
  )

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
# expect_that(pData(result),is_a("data.table"))
expect_that(colnames(pData(result)),equals(c("UID","TCELLSUBSET","CYTOKINE","RefTreat","ANTIGEN")))
expect_that(nrow(pData(result)),equals(150))


