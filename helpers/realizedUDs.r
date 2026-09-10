# To run numerical simulations  of 1-d switching Langevin:
# automated, based on master.r

# to convert each habitat covariate into functions that directly give
# values of covariate and of gradient (model-specific)
num_fns <- preprocess(model$K,controls,habitat)

# calculate UDs implied by within-state selection and marginal state
# to allow initialisation and plotting
UD_out <- UDs(model,controls)
if (is.null(model$beta_zero)) model$beta_zero <- UD_out$derived_beta_zero

# carry out simulations

# naive version
model$Qfixed <- TRUE
if (!is.null(controls$seed)) set.seed(controls$seed)
outputs_naive <- iterate(model,controls,UD_out,num_fns)

# corrected version
model$Qfixed <- FALSE
if (!is.null(controls$seed)) set.seed(controls$seed)
outputs_novel <- iterate(model,controls,UD_out,num_fns)

save(UD_out,outputs_naive,outputs_novel,
     file=paste0(scenario,"/outputs_beta",paste0(betaCoeff,collapse="_"),"_sigma",paste0(sig,collapse="_"),".RData"))

source("helpers/plotRealizedUDs.r")
