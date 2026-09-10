# Simulate given model, controls
simulate <- function(model, controls,UD_out,num_fns)
{
  S <- model$S
  K <- model$K
  sigma <- model$sigma
  beta <- model$beta

  dt <- controls$dt
  tt <- seq(0,controls$max_t,dt)
  nt <- length(tt)

  # Drift leading term
  v <- (model$sigma^2/2)*dt

  # Noise multiplier
  m <- model$sigma*sqrt(dt)

  # To save over time
  ss <- xx <- rep(NA,nt)

  # To overwrite along the way
  g <- matrix(NA,nt,K)

  eps <- rnorm(nt)

  initial_values <- initialize_sim(model,controls,UD_out$joint_ud,UD_out$conditional_ud,n=NULL) # n is NULL if N=1; in fact, for now, always NULL!
  ss[1] <- s <- initial_values$s
  xx[1] <- x <- initial_values$x

  t <- 0
  tt <- (0:(nt-1))*dt

  Qfixed <- model$Qfixed
  if (Qfixed)
  {
    # Transition probabilities
    P <- Matrix::expm(model$Q*dt)
    # Pre-simulate MC!
    u <- runif(nt)
    PP <- apply(cbind(0,P),1,cumsum)
    poss <- .bincode(u,PP[,1])
    for (s in 2:S) poss <- cbind(poss,.bincode(u,PP[,s]))
  }

  h <- matrix(NA,nt,K)

  # Main loop: i=1 for starting values
  for(i in 2:nt)
  {
    for (k in 1:K)
      g[i,k] <- num_fns$grad_num[[k]](x)
    drift <- v[s]*sum(beta[s,]*g[i,])
    noise <- m[s]*eps[i]
    new_x <- x + drift + noise
    if (Qfixed) {
      new_s <- poss[i,s]
      for (k in 1:K)
        h[i,k] <- num_fns$habi_num[[k]](x)
    } else
    {
      for (k in 1:K)
        h[i,k] <- num_fns$habi_num[[k]](x)

      eta <- as.vector(model$beta_zero + model$beta %*% h[i,])

      Qmat <- matrix(0, S, S)
      for (state in 1:S) {
        Qmat[state, ] <- model$rate_function(eta, state, model)
        Qmat[state, state] <- 0
        Qmat[state, state] <- -sum(Qmat[state, ])
      }

      if (S == 2) {
        lambda <- Qmat[1, 2] + Qmat[2, 1]
        exp_lt <- exp(-lambda * dt)
        P_row <- numeric(2)
        if (lambda > 0) {
          if (s == 1) {
            P_row <- c((Qmat[2, 1] + Qmat[1, 2] * exp_lt) / lambda, (Qmat[1, 2] - Qmat[1, 2] * exp_lt) / lambda)
          } else {
            P_row <- c((Qmat[2, 1] - Qmat[2, 1] * exp_lt) / lambda, (Qmat[1, 2] + Qmat[2, 1] * exp_lt) / lambda)
          }
        } else {
          P_row[s] <- 1
        }
        p <- P_row
      } else {
        # Fallback to Matrix::expm for generalized S > 2 models
        Pmat <- as.matrix(Matrix::expm(Qmat * dt))
        p <- Pmat[s, ]
      }

      new_s <- sample(1:S,prob=p,size=1)
    }
    xx[i] <- x <- new_x
    ss[i] <- s <- new_s
  }
  return(list(xx=xx,ss=ss,tt=tt,g=g,h=h))
}
