# Set model parameters

S <- 2 # number of states
K <- 2 # number of covariates

q1 <- 3
q2 <- 5
rho <- 0.01

model <- list(
  S=S,
  K=K,
  Q=matrix(c(-q1,q1,q2,-q2)*rho,S,S,byrow=TRUE),
  Qfixed=TRUE,
  Psi=matrix(c(0,1,1,0)*sqrt(q1*q2)*rho,S,S),
  #rate_function=MetropolisHastings_rates,
  #rate_function=Gibbs_rates,
  #rate_function=NonCanonical_rates,
  rate_function=SymmetricRoot_rates,
  sigma=sig,
  beta_zero=NULL, #c(log(5/3),0),
  beta=matrix(betaCoeff,
              S,K,byrow=TRUE) # state then covariate
)

#model$header <- if(model$Qfixed) "Fixed transition rates" else "Targetted transition rates"
rm(S,K,q1,q2,rho)
