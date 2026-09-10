initialize_sim <- function(model,controls,joint_ud,conditional_ud,n=1) # n is NULL if N=1
{
  s <- controls$s0
  x <- controls$x0
  if (is.null(x))
  {
    # No initial x value given
    if (is.null(s))
    {
      # Construct marginal for s, then sample
      marginal_state <- rowSums(joint_ud)
      #cat("Marginal for state",marginal_state/sum(marginal_state),"\n")
      s <- sample(1:model$S,size=1,prob=marginal_state)
    } 
    # Using given or sampled s,construct x
    xx <- seq(controls$xmin,controls$xmax,controls$dx)
    x <- sample(xx,size=1,prob=conditional_ud[s,])+runif(1,-controls$dx/2,controls$dx/2)
  } else # x given
    if (is.null(s))
    {
      # Only initial x given
      x_index <- round((x-controls$xmin)/controls$dx)+1
      conditional_state <- joint_ud[,x_index]
      s <- sample(1:model$S,size=1,prob=conditional_state)
    } # else both initial values given, so nothing to do
  return(list(s=s,x=x))
}