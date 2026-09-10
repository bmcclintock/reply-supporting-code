initial_plots <- function(model,controls,UD_out)
{
xx <- seq(controls$xmin,controls$xmax,controls$dx)
nx <- length(xx)

S <- model$S
K <- model$K

# Optionally, plot covariates
#hrange <- range(c(hh))  
#if (FALSE)
#{
#plot(xx,rep(0,nx),ylim=hrange,
#     xlab="x",
#     ylab="Covariate",     
#     type="l",lty=3)
#for(k in 1:K)
#  lines(xx,hh[k,]) 
#}


umax <- max(c(UD_out$conditional_ud))
#pdf("UD.pdf",width=10)

magic <- if (is.null(appear$magic)) 1.2 else appear$magic
# Simply a "magic number" that replaces the 1.05 default in R, to give some wiggle room in plot
plot(xx,rep(0,nx),
     xlim=c(appear$xmin,appear$xmax),xaxs="i",
     ylim=c(0,umax*magic),yaxs="i",
     xlab="x",
     ylab="Utilization Distributions",
     #main=model$header,
     type="l",lty=3)
for (s in 1:S)
{
  lines(xx,UD_out$conditional_ud[s,],
        col=appear$state_colours[s],lwd=2,lty=5)
}
lines(xx,UD_out$marginal_ud,col=1,lwd=2,lty=5)
return(invisible(NULL)) # exists only for side-effect of plotting
}