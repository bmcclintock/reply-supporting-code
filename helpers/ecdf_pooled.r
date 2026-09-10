ecdf_pooled <- function(model,controls,appear,outputs,out_pooled,UD_out,header)
{
S <- model$S

# The x values used to get numerical UDs
xx <- seq(controls$xmin,controls$xmax,controls$dx)

plot(xx,cumsum(UD_out$marginal_ud)*controls$dx,
     xlab="x",ylab="Cumulative probability",
     main=header,
     lty=5,lwd=2,type="l",xlim=c(appear$xmin,appear$xmax))
for (s in 1:S)
{
  lines(xx,cumsum(UD_out$conditional_ud[s,])*controls$dx,lty=5,
        col=appear$state_colours[s],lwd=2)
}

thin <- 50

# Each run separately

if (appear$show_runs_ecdf)
{
indices <- seq(1,length(outputs[[1]]$xx),by=thin)
for (n in 1:controls$N)
{
  x_thinned <- outputs[[n]]$xx[indices]
  s_thinned <- outputs[[n]]$ss[indices]
  plot(ecdf(x_thinned),add=TRUE,
       do.points=FALSE,verticals=TRUE,lwd=1)
  for (s in 1:S)
  {
    plot(ecdf(x_thinned[s_thinned==s]),add=TRUE,
         do.points=FALSE,verticals=TRUE,
         col=appear$state_colours[s],lwd=1)
  }
}
}

# Pooled runs - assuming display_pooled has already been run
# to create xx_pooled, ss_pooled

indices <- seq(1,length(out_pooled$xx_pooled),by=thin)
xx_thinned <- out_pooled$xx_pooled[indices]
ss_thinned <- out_pooled$ss_pooled[indices]
plot(ecdf(xx_thinned),add=TRUE,
     do.points=FALSE,verticals=TRUE,lwd=2)
for (s in 1:S)
{
  plot(ecdf(xx_thinned[ss_thinned==s]),add=TRUE,
       do.points=FALSE,verticals=TRUE,
       col=appear$state_colours[s],lwd=2)
}
return(invisible(NULL))
}