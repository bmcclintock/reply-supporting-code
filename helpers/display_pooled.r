display_pooled <- function(S,controls,appear,outputs,header)
{
xx_pooled <- ss_pooled <- c()
N <- controls$N
for (n in 1:N)
{
  xx_pooled <- c(xx_pooled,outputs[[n]]$xx)
  ss_pooled <- c(ss_pooled,outputs[[n]]$ss)
  usage_overall <- density(outputs[[n]]$xx,
                           bw=appear$xbw,from=appear$xmin,to=appear$xmax,kernel="triangular")
  usage_by_state <- list()
  for (s in 1:S)
  {
    usage_by_state[[s]] <- density(outputs[[n]]$xx[outputs[[n]]$ss==s],
                                   bw=appear$xbw,from=appear$xmin,to=appear$xmax,kernel="triangular")
  }
  
  if (appear$show_runs_kde)
  {
  lines(usage_overall,lwd=1)
  for (s in 1:S)
    lines(usage_by_state[[s]],lwd=1,col=appear$state_colours[s])
  }
}

usage_pooled <- density(xx_pooled,
                         bw=appear$xbw,from=appear$xmin,to=appear$xmax,kernel="triangular")

usage_pooled_by_state <- list()
for (s in 1:S)
{
  usage_pooled_by_state[[s]] <- density(xx_pooled[ss_pooled==s],
                                        bw=appear$xbw,from=appear$xmin,to=appear$xmax,kernel="triangular")
}

lines(usage_pooled,lwd=2)
#title(paste0(model$header,", KDE bandwidth ",appear$xbw))#,line=1)
title(paste0(header,", KDE bandwidth ",appear$xbw))#,line=1)
for (s in 1:S)
  lines(usage_pooled_by_state[[s]],lwd=2,col=appear$state_colours[s])
  
return(list(xx_pooled=xx_pooled,ss_pooled=ss_pooled))
}
