# Define environment for each k=1,...,K
habitat <- list(
  list(fun="laplace_patch",left=-15,right=0,sigma=grad_scale),
  list(fun="laplace_patch",left=5,right=15,sigma=grad_scale)
)
