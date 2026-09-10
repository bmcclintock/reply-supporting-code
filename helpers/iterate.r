iterate <- function(model,controls,UD_out,num_fns)
{
  N <- controls$N
  outputs <- list()
  for (n in 1:N)
  {
    outputs[[n]]<-simulate(model,controls,UD_out,num_fns)
  }
  return(outputs)
}