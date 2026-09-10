Gibbs_rates <- function(eta,s,model)
{
  rho <- model$rho
  target <- exp(eta)
  return(rho*target/sum(target))
}

MetropolisHastings_rates <- function(eta,s,model)
{
  Q <- model$Q
  return(Q[s,]*pmin(1,exp(eta-eta[s])*Q[,s]/Q[s,]))
}

NonCanonical_rates <- function(eta,s,model)
{
  Q <- model$Q
  target <- exp(eta)
  numerator <- target*Q[,s]*Q[s,]
  denominator <- target[s]*Q[s,]+target*Q[,s]
  return(numerator/denominator)
}

SymmetricRoot_rates <- function(eta,s,model)
{
  if (is.null(model$Psi)) Psi <- sqrt(model$Q*t(model$Q)) else Psi <- model$Psi
  target <- exp(eta)
  return(Psi[s,]*sqrt(target/target[s]))
}