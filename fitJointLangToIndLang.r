source("helpers/jointLangevin.R")

scenario <- "TwoPatchRepel"
if(!dir.exists(scenario)) dir.create(scenario)

nSims <- 100

betaCoeff <- c(2,-0.5,-0.5,2) # habitat selection coefficients
sig <- c(1,1) # speed parameter
beta0 <- matrix(c(log(0.03),log(0.05)),1,2) # state transition rates (working scale)
grad_scale <- 5 # gradient scaling factor

statesKnown <- FALSE

source("helpers/setup.r")
source(paste0(scenario,"/model.r"))
source(paste0(scenario,"/habitat.r"))
source(paste0(scenario,"/controls.r"))
source(paste0(scenario,"/appear.r"))

seed <- 0

# Modified to accommodate CIs derived from the Hessian natively
xResults_est <- matrix(NA, nSims, 4)
xResults_lower <- matrix(NA, nSims, 4)
xResults_upper <- matrix(NA, nSims, 4)

sigResults_est <- matrix(NA, nSims, 2)
sigResults_lower <- matrix(NA, nSims, 2)
sigResults_upper <- matrix(NA, nSims, 2)

stateResults <- timeInState1 <- rep(NA,nSims)

for(isim in 1:nSims){
  outputs_novel <- st <- tryCatch(stop(),error=function(e) e)

  message("Simulation ",isim)
  while(inherits(outputs_novel,"error") | inherits(st,"error")){
    source("helpers/setup.r")
    source(paste0(scenario,"/controls.r"))
    source(paste0(scenario,"/model.r"))
    source(paste0(scenario,"/habitat.r"))

    controls$seed <- seed
    num_fns <- preprocess(model$K,controls,habitat)
    UD_out <- UDs(model,controls)
    if (is.null(model$beta_zero)) model$beta_zero <- UD_out$derived_beta_zero

    # simulate from the independent model
    model$Qfixed <- TRUE
    if (!is.null(controls$seed)) set.seed(controls$seed)
    outputs_novel <- tryCatch(iterate(model,controls,UD_out,num_fns),error=function(e) e)

    if(!inherits(outputs_novel,"error")){

      x_obs <- unlist(lapply(outputs_novel, function(x) x$xx))
      g1_obs <- unlist(lapply(outputs_novel, function(x) x$g[,1]))
      g2_obs <- unlist(lapply(outputs_novel, function(x) x$g[,2]))
      h1_obs <- unlist(lapply(outputs_novel, function(x) x$h[,1]))
      h2_obs <- unlist(lapply(outputs_novel, function(x) x$h[,2]))
      states_true <- unlist(lapply(outputs_novel, function(x) x$ss))

      track_start <- as.integer(is.na(g1_obs))

      g1_safe <- ifelse(is.na(g1_obs), 0, g1_obs)
      g2_safe <- ifelse(is.na(g2_obs), 0, g2_obs)
      h1_safe <- ifelse(is.na(h1_obs), 0, h1_obs)
      h2_safe <- ifelse(is.na(h2_obs), 0, h2_obs)

      par <- c(model$beta_zero,model$beta,log(model$sigma),log(model$Psi[2]),0)

      fit <- nlminb(par, fitjoint_cpp,
                                x=x_obs, g1=g1_safe, g2=g2_safe,
                                h1=h1_safe, h2=h2_safe,
                                dt=controls$dt, track_start=track_start)

      hess <- optimHess(fit$par, fitjoint_cpp,
                        x=x_obs, g1=g1_safe, g2=g2_safe,
                        h1=h1_safe, h2=h2_safe,
                        dt=controls$dt, track_start=track_start)

      se <- tryCatch(sqrt(diag(solve(hess))), error = function(e) rep(NA, 10))

      st <- tryCatch(jointviterbi_cpp(fit$par, x_obs, g1_safe, g2_safe,
                                      h1_safe, h2_safe, controls$dt, track_start),
                     error=function(e) e)

      # check for label switching
      if (mean(st == states_true) < 0.5) {

        fitpar <- fit$par
        fitse <- se

        fit$par[1:2] <- fitpar[c(2, 1)]
        se[1:2] <- fitse[c(2, 1)]

        fit$par[3:4] <- fitpar[c(4, 3)]
        se[3:4] <- fitse[c(4, 3)]

        fit$par[5:6] <- fitpar[c(6, 5)]
        se[5:6] <- fitse[c(6, 5)]

        fit$par[7:8] <- fitpar[c(8, 7)]
        se[7:8] <- fitse[c(8, 7)]

        if (!inherits(st, "error")) {
          st <- 3 - st
        }
      }

      xResults_est[isim,] <- fit$par[3:6]
      xResults_lower[isim,] <- fit$par[3:6] - 1.96 * se[3:6]
      xResults_upper[isim,] <- fit$par[3:6] + 1.96 * se[3:6]

      sigResults_est[isim,] <- exp(fit$par[7:8])
      sigResults_lower[isim,] <- exp(fit$par[7:8] - 1.96 * se[7:8])
      sigResults_upper[isim,] <- exp(fit$par[7:8] + 1.96 * se[7:8])

      if(!inherits(st,"error")){
        timeInState1[isim] <- sum(st == 1) / length(st)
        stateResults[isim] <- mean(st == states_true)
      }
    }
    seed <- seed+1
  }
  print(c(apply(xResults_est,2,mean,na.rm=TRUE),apply(sigResults_est,2,mean,na.rm=TRUE),mean(stateResults,na.rm=TRUE)))
  vp <- cbind(xResults_est, sigResults_est, stateResults)
  colnames(vp) <- c("beta_1,0","beta_1,1","beta_2,0","beta_2,1","sigma_1","sigma_2","viterbi")
  vioplot::vioplot(vp,ylim=c(-0.5,3))
  truevals <- c(betaCoeff,sig,0.95)
  for(i in 1:length(truevals)){
    lines(c(i-1,i)+0.5,rep(truevals[i],2),col="red",lty=2,lwd=2)
  }
}

summary(stateResults)

# Fixed extraction logic for matrix structures adding coverage probability column
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
pdf(paste0(scenario,"/JointToIndResults_beta",paste0(betaCoeff,collapse="_"),"_sigma",paste0(sig,collapse="_"),'_gradscale',grad_scale,ifelse(statesKnown,"_statesKnown",""),".pdf"),width=11,height=8)
vioplot::vioplot(vp,ylim=c(-0.5,3))
truevals <- c(betaCoeff,sig)
for(i in 1:length(truevals)){
  lines(c(i-1,i)+0.5,rep(truevals[i],2),col="red",lty=2,lwd=2)
}
dev.off()

save(stateResults, timeInState1, xdf, xResults_est, sigResults_est, truevals, grad_scale,
     file=paste0(scenario,"/JointToInd_beta",paste0(betaCoeff,collapse="_"),"_sigma",paste0(sig,collapse="_"),'_gradscale',grad_scale,ifelse(statesKnown,"_statesKnown",""),".RData"))
