# Calculate true distributions to allow sampling of initial state and location
UDs <- function(model,controls)
{
K <- model$K

# Set up x values to consider
xx <- seq(controls$xmin,controls$xmax,controls$dx)
nx <- length(xx)

# Evaluate covariates
hh <- matrix(NA,K,nx)
for(k in 1:K)
  hh[k,] <- eval(call(habitat[[k]]$fun,xx,habitat[[k]]))

# Calculate the long-run UDs for each individual state
# Use same xx and hh as above
S <- model$S
conditional_ud <- joint_ud <- matrix(NA,S,nx)

if (is.null(model$beta_zero))
{
  # beta_zero should be set based on Q
  if (model$S==2) 
  {
    # easy to get stationary probabilities for naive approach
    Q <- model$Q
    p <- c(Q[2,2],Q[1,1])/((Q[1,1]+Q[2,2]))
    # calculate conditionals first
    sum_eLP <- rep(NA,S)
    for (s in 1:S)
    {
      # LP is linear predictor, or eta excluding beta_zero
      LP <- rep(0,nx)
      for (k in 1:K)
      {
        LP <- LP + model$beta[s,k]*hh[k,]
      }
      eLP <- exp(LP)
      sum_eLP[s] <- sum(eLP)
      conditional_ud[s,] <- eLP/sum_eLP[s]/controls$dx
      joint_ud[s,] <- conditional_ud[s,]*p[s]
    }
    # Equivalently
    #model$beta_zero <- log(p/sum_eLP)
    derived_beta_zero <- log(p/sum_eLP)
  } else
    stop("Need stationary distribution for S>2")
} else
{
  # beta_zero given
  for (s in 1:S)
  {
    # LP is linear predictor, or eta
    LP <- rep(model$beta_zero[s],nx)
    for (k in 1:K)
    {
      LP <- LP + model$beta[s,k]*hh[k,]
    }
    eLP <- exp(LP)
    joint_ud[s,] <- eLP
    conditional_ud[s,] <- eLP/sum(eLP)/controls$dx
  }
}

# Correct marginal UD
marginal_ud <- colSums(joint_ud)/sum(joint_ud)/controls$dx

#rm(eLP,k,K,LP,nx,p,Q,s,S,sum_eLP,xx)

# modifies model$beta_zero (if NULL) 
# defines conditional_ud,joint_ud,marginal_ud,hh
return(list(conditional_ud=conditional_ud,joint_ud=joint_ud,marginal_ud=marginal_ud,hh=hh,derived_beta_zero=derived_beta_zero))
}