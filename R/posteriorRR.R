#'Compute the posterior response rate from MCMC samples
#'
#'@description 
#'Uses the sampled response indicator when the model class is MCMCResult (or MIMOSAResult encapsulating an MCMCResult) to
#'compute the posterior response rate distribution. This is summarized to its median, 2.5th and 97.5th credible interval. 
#'
#'@details 
#'The posterior response rate is the correct way to compare response rates across studies and treatment groups, as it is unbiased compared to the response rate computed from hard thresholds
#'of posterior response probabilities. The credible intervals have the correct behaviour as sample size increases. 
#'
#'Future versions will allow passing of quantiles for summarization, and maybe the full distribution, depending on needs.
#'@param x A MIMOSAResultList or MIMOSAResult. All models should be of type MCMCResult, or fitted using method="mcmc" in MIMOSA.
#'@param variable an unquoted variable name in the pData() table of all the models that specifies the grouping variable for which to compute response rates.
#'@return a tibble with the grouping variable, and the 2.5th, median, and 97.5th percentiles of the posterior response rate.
#'@import tidyr
#'@importFrom dplyr bind_rows bind_cols group_by summarize ungroup `%>%`
# @importFrom plyr ldply
#'@export
getPosteriorResponseRate <- function(x, variable) {
  quo_variable <- enquo(variable)
  if ("MIMOSAResultList" %in% class(x)) {
      retme<-lapply(x,function(x)getPosteriorResponseRate(x,!!quo_variable))
      retme<-plyr::ldply(retme,.id="model")
#      retme2<-dplyr::bind_rows(retme)
#      retme<-dplyr::bind_cols(retme2,model=names(retme))
      retme
  } else if ("MIMOSAResult" %in% class(x)& "MCMCResult"%in%class(x@result)) {
    mat <- x@result@IndMat
    cn <- data.frame(rn = colnames(mat))
    pd <- pData(x)
    groups <-
      pd %>% dplyr::select(!!quo_variable) %>% dplyr::bind_cols(rn = rownames(pd))
    posterior_rr<-dplyr::right_join(groups, cn) %>% dplyr::group_by(!!quo_variable) %>%
      dplyr::do(., {
        data.frame(rr = rowMeans(mat[, .$rn])) %>% dplyr::summarize(
          Q2.5 = quantile(rr, 0.025),
          Q50 = quantile(rr, 0.5),
          Q97.5 = quantile(rr, 0.975)
        )
      })
    posterior_rr%>%dplyr::ungroup()
  }
}
