# Fit
fit <- function(model,habitat,bin_out)
{
# In bin_out, have: 
#bin_marginal_obs[x]
#bin_conditional_obs[s,x]
#x_cuts
#n_obs
#n_obs_by_state[s]

K <- model$K
S <- model$S

# Admin
x_cuts <- bin_out$x_cuts
nc <- length(x_cuts)
bin_x <- (x_cuts[-1]+x_cuts[-nc])/2

# Temporary fudge to avoid having to "construct" formula in GLM
stopifnot(K==2) 

# Evaluate covariates
bin_hh <- matrix(NA,K,nc-1)
for(k in 1:K)
  bin_hh[k,] <- eval(call(habitat[[k]]$fun,bin_x,habitat[[k]]))

# Marginal
bin_marginal_counts <- bin_out$bin_marginal_obs*bin_out$n_obs

glm0 <- glm(bin_marginal_counts~bin_hh[1,]+bin_hh[2,],family="quasipoisson")

# Conditional
bin_conditional_counts <- matrix(NA,model$S,nc-1)
glm_s <- list()

for (s in 1:S)
{
bin_conditional_counts[s,] <- bin_out$bin_conditional_obs[s,]*bin_out$n_obs_by_state[s]
glm_s[[s]] <- glm(bin_conditional_counts[s,]~bin_hh[1,]+bin_hh[2,],family="quasipoisson")
}

return(list(glm0=glm0,glm_s=glm_s))
}

