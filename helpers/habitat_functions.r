# Generic function to get gradient at a single point
gradient1 <- function(x,k,habitat,delta=10^-3)
{
  hh <- eval(call(habitat[[k]]$fun,x+c(-1,1)*delta,habitat[[k]]))
  return(diff(hh)/(2*delta))
}

# Functions to define covariates

# Simple functions
zero <- function(x,...) 
  return(rep(0,length(x)))

sinusoid <- function(x,settings) 
  return(sin((x/settings$period+settings$phase)*2*pi)*settings$amplitude)

quadratic <- function(x,settings)
  return(settings$a*x^2+settings$b*x+settings$c)

# Kernels of distributions, on log scale, with max=0

normal <- function(x,settings)
  return(-(x-settings$mu)^2/(2*settings$sigma^2))

laplace <- function(x,settings)
  return(-abs(x-settings$mu)/(2*settings$sigma))

logistic <- function(x,settings)
  return(dlogis(x,location=settings$mu,scale=settings$sigma,log=TRUE)-dlogis(0,0,scale=settings$sigma,log=TRUE))

# Multiple "point" patches

patches_linear <- function(x,settings)
{
  d <- abs(outer(x,settings$nuclei,"-"))
  return(-apply(d,1,min))
}

# Multiple "point" patches with kernels as above
# Vector of nuclei replaces mu
# Scale sigma assumed common across patches

normal_nuclei <- function(x,settings)
{
  d <- outer(x,settings$nuclei,"-")
  return(-apply(d^2,1,min)/(2*settings$sigma^2))
}

laplace_nuclei <- function(x,settings)
  # Scaled version of patches_linear
{
  d <- abs(outer(x,settings$nuclei,"-"))
  return(-apply(d,1,min)/settings$sigma)
}

logistic_nuclei <- function(x,settings)
  # Scaled version of patches_linear
{
  #print(summary(x))
  #print(settings$nuclei)
  d <- abs(outer(x,settings$nuclei,"-"))
  #print(summary(c(d)))
  #print(summary(apply(d,1,min)))
  #print(settings$sigma)
  return(dlogis(apply(d,1,min),location=0,scale=settings$sigma,log=TRUE)-dlogis(0,0,scale=settings$sigma,log=TRUE))
}

laplace_patch <- function(x,settings)
{
  return(-(pmax(x-settings$right,0)+pmax(settings$left-x,0))/settings$sigma)
}

laplace_patches <- function(x,settings)
{
  dR <- pmax(outer(x,settings$right,"-"),0)
  dL <- -pmin(outer(x,settings$left,"-"),0)
  d <- dR+dL
  return(-apply(d,1,min)/settings$sigma)
}
