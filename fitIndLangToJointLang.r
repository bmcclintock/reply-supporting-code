source("helpers/indLangevin.R")

scenario <- "TwoPatchRepel"
if(!dir.exists(scenario)) dir.create(scenario)

nSims <- 100

betaCoeff <- c(2,-0.5,-0.5,2) # habitat selection coefficients
sig <- c(1,1) # speed parameter
beta0 <- matrix(-5,1,2) # initial values for state transition rates
grad_scale <- 5 # gradient scaling factor

statesKnown <- FALSE

source("helpers/setup.r")
source(paste0(scenario,"/model.r"))
source(paste0(scenario,"/habitat.r"))
source(paste0(scenario,"/controls.r"))
source(paste0(scenario,"/appear.r"))

GLMs <- TRUE
PDFs <- TRUE
destination <- paste0(scenario,"/figures/")
if(!dir.exists(destination)) dir.create(destination)
appear$xmin <- -35
appear$xmax <- 30
source("helpers/realizedUDs.r")

seed <- 0

xResults_est <- matrix(NA, nSims, 4)
xResults_lower <- matrix(NA, nSims, 4)
xResults_upper <- matrix(NA, nSims, 4)

sigResults_est <- matrix(NA, nSims, 2)
sigResults_lower <- matrix(NA, nSims, 2)
sigResults_upper <- matrix(NA, nSims, 2)

stateResults <- timeInState1 <- rep(NA, nSims)

# Initial parameters for optimization (beta1_1, beta1_2, beta2_1, beta2_2, log_sig1, log_sig2, log_q12, log_q21)
init_par <- c(betaCoeff, log(sig), beta0)

for(isim in 1:nSims){
  outputs_novel <- tryCatch(stop(),error=function(e) e)

  message("Simulation ",isim)
  while(inherits(outputs_novel,"error")){
    source("helpers/setup.r")
    source(paste0(scenario,"/controls.r"))
    source(paste0(scenario,"/model.r"))
    source(paste0(scenario,"/habitat.r"))

    controls$seed <- seed
    num_fns <- preprocess(model$K,controls,habitat)

    UD_out <- UDs(model,controls)
    if (is.null(model$beta_zero)) model$beta_zero <- UD_out$derived_beta_zero

    # simulate from the joint model
    model$Qfixed <- FALSE
    if (!is.null(controls$seed)) set.seed(controls$seed)
    outputs_novel <- tryCatch(iterate(model,controls,UD_out,num_fns),error=function(e) e)

    if(!inherits(outputs_novel,"error")){

      x_obs <- unlist(lapply(outputs_novel, function(x) x$xx))
      g1_obs <- unlist(lapply(outputs_novel, function(x) x$g[,1]))
      g2_obs <- unlist(lapply(outputs_novel, function(x) x$g[,2]))
      states_true <- unlist(lapply(outputs_novel, function(x) x$ss))

      track_start <- as.integer(is.na(g1_obs))

      g1_safe <- ifelse(is.na(g1_obs), 0, g1_obs)
      g2_safe <- ifelse(is.na(g2_obs), 0, g2_obs)

      fit <- optim(init_par, langevin_nll_cpp,
                   x = x_obs, g1 = g1_safe, g2 = g2_safe,
                   dt = controls$dt, track_start = track_start,
                   method = "BFGS", hessian = TRUE)

      se <- tryCatch(sqrt(diag(solve(fit$hessian))), error = function(e) rep(NA, 8))

      xResults_est[isim,] <- fit$par[1:4]
      xResults_lower[isim,] <- fit$par[1:4] - 1.96 * se[1:4]
      xResults_upper[isim,] <- fit$par[1:4] + 1.96 * se[1:4]

      sigResults_est[isim,] <- exp(fit$par[5:6])
      sigResults_lower[isim,] <- exp(fit$par[5:6] - 1.96 * se[5:6])
      sigResults_upper[isim,] <- exp(fit$par[5:6] + 1.96 * se[5:6])

      st <- viterbi_cpp(fit$par, x_obs, g1_safe, g2_safe, controls$dt, track_start)

      timeInState1[isim] <- sum(st == 1) / length(st)
      stateResults[isim] <- mean(st == states_true)
    }
    seed <- seed + 1
  }
  print(c(apply(xResults_est, 2, mean, na.rm = TRUE), apply(sigResults_est, 2, mean, na.rm = TRUE), mean(stateResults, na.rm = TRUE)))
}

summary(stateResults)

xdf <- cbind(t(apply(xResults_est, 2, summary)),
             colMeans(xResults_lower <= matrix(betaCoeff, nSims, 4, byrow=TRUE) &
                        xResults_upper >= matrix(betaCoeff, nSims, 4, byrow=TRUE), na.rm=TRUE))

sig_df <- cbind(t(apply(sigResults_est, 2, summary)),
                colMeans(sigResults_lower <= matrix(sig, nSims, 2, byrow=TRUE) &
                           sigResults_upper >= matrix(sig, nSims, 2, byrow=TRUE), na.rm=TRUE))

xdf <- rbind(xdf, sig_df)
colnames(xdf)[7] <- "Coverage"

vp <- cbind(xResults_est, sigResults_est)
colnames(vp) <- c("beta_1,1","beta_1,2","beta_2,1","beta_2,2","sigma_1","sigma_2")

pdf(paste0(scenario,"/IndToJointResults_beta",paste0(betaCoeff,collapse="_"),"_sigma",paste0(sig,collapse="_"),'_gradscale',grad_scale,ifelse(statesKnown,"_statesKnown",""),".pdf"),width=11,height=8)
vioplot::vioplot(vp)
truevals <- c(betaCoeff,sig)
for(i in 1:length(truevals)){
  lines(c(i-1,i)+0.5,rep(truevals[i],2),col="red",lty=2,lwd=2)
}
dev.off()

save(stateResults, timeInState1, xdf, xResults_est, sigResults_est, truevals, grad_scale,
     file=paste0(scenario,"/IndToJoint_beta",paste0(betaCoeff,collapse="_"),"_sigma",paste0(sig,collapse="_"),'_gradscale',grad_scale,ifelse(statesKnown,"_statesKnown",""),".RData"))
