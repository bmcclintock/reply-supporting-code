# to convert each habitat covariate into functions that directly give 
# values of covariate and of gradient (model-specific)
preprocess <- function(K,controls,habitat)
{
delta <- if(is.null(controls$delta)) 10^-3 else controls$delta
K <- model$K
xx <- seq(controls$xmin,controls$xmax,by=controls$dx)
habi_num <- grad_num <- list()

for (k in 1:K)
{
  hh <- eval(call(habitat[[k]]$fun,xx,habitat[[k]]))
  habi_num[[k]] <- approxfun(xx,hh,rule=2)
  h_minus <- eval(call(habitat[[k]]$fun,xx-delta,habitat[[k]]))
  h_plus <- eval(call(habitat[[k]]$fun,xx+delta,habitat[[k]]))
  gg <- (h_plus-h_minus)/(2*delta)
  grad_num[[k]] <- approxfun(xx,gg,rule=2)
}
#rm(delta,gg,hh,h_minus,h_plus,k,K,xx)
# Creates habi_num,grad_num
return(list(habi_num=habi_num,grad_num=grad_num))
}