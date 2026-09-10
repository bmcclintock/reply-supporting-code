# Probs in bins
bins <- function(S,controls,appear,UD_out,output_list,header)
{
# stopifnot(controls$N>1) ## should work anyway

xx_out <- output_list$xx_pooled
ss_out <- output_list$ss_pooled


# The x values used to get numerical UDs
xx <- seq(controls$xmin,controls$xmax,controls$dx)
nx <- length(xx)

# Settings for the bins
bin_width <- if (is.null(appear$bin_width)) 1.0 else appear$bin_width
x_cuts <- bin_width*seq(floor(min(xx_out)/bin_width),ceiling(max(xx_out)/bin_width),by=1)
nbin <- length(x_cuts)-1

# The marginal case
# CDF
marginal_CDF <- approxfun(y=cumsum(UD_out$marginal_ud)*controls$dx,x=xx,
                          yleft=0,yright=1,rule=2)
# Bin probabilities
n_obs <- length(xx_out)
bin_marginal_probs <- diff(marginal_CDF(x_cuts))
bin_marginal_obs <- hist(xx_out,breaks=x_cuts,plot=FALSE)$counts/n_obs

# The conditional case, for each state
bin_conditional_probs <- bin_conditional_obs <- matrix(NA,S,nbin)
n_obs_by_state <- rep(NA,S)
for (s in 1:S)
{
  n_obs_by_state[s] <- sum(ss_out==s)
  conditional_CDF <- approxfun(y=cumsum(UD_out$conditional_ud[s,])*controls$dx,x=xx,
                               yleft=0,yright=1,rule=2)
  bin_conditional_probs[s,] <- diff(conditional_CDF(x_cuts))
  bin_conditional_obs[s,] <- hist(xx_out[ss_out==s],breaks=x_cuts,plot=FALSE)$counts/n_obs_by_state[s]
}

# Set limits to accommodate all cases
bin_max  <- if(is.null(appear$bin_max)) 0 else appear$bin_max
limits <- range(bin_conditional_obs,bin_conditional_probs,
                bin_marginal_obs,bin_marginal_probs,0,bin_max)
# Prettyfying...
limits <- range(pretty(limits))

# Captions
x_cap <- "Theoretical probabilities"
y_cap <- "Empirical probabilities"


# Plot marginal case
plot(bin_marginal_probs,bin_marginal_obs,pch=3,asp=1,xlim=limits,ylim=limits,
     xaxs="i",yaxs="i",
     main=paste("Cell probabilities, cell size",bin_width,"\n",model$header),
     xlab=x_cap,ylab=y_cap)
abline(0,1,lty=3)

# Plot conditional case for each state
for (s in 1:S)
{
  plot(bin_conditional_probs[s,],bin_conditional_obs[s,],
       pch=3,asp=1,xlim=limits,ylim=limits,
       xaxs="i",yaxs="i",
       xlab=x_cap,ylab=y_cap,
       #main=paste(model$header,"\n","State",s),
       main=paste("\nState",s),
       col=appear$state_colours[s])
  abline(0,1,lty=3)
}
return(list(x_cuts=x_cuts,bin_marginal_obs=bin_marginal_obs,bin_conditional_obs=bin_conditional_obs,n_obs=n_obs,n_obs_by_state=n_obs_by_state))
}
