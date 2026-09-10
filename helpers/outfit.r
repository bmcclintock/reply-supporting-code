outfit <- function(model,glms,outfile,scenario,header)
{
# Assuming K=2!
stopifnot(model$K==2)

cat(scenario,"\n\n",
    file=outfile,append=TRUE)

# Overall
summ <- summary(glms$glm0,correlation=TRUE)
cat(header,"\n",
    "Fitted:\n",
    "Beta 1:",summ$coeff[2,1:2],"\n",
    "Beta 2:",summ$coeff[3,1:2],"\n",
    "Dispersion",summ$dispersion,"\n",
    "Degrees of freedom",summ$df[2],"\n",
    "Correlation",summ$correlation[3,2],"\n\n",
    file=outfile,append=TRUE)

#conditional
for (s in 1:model$S)
{
  summ <- summary(glms$glm_s[[s]],correlation=TRUE)
  cat("s =",s,"\n",
      "True beta 1:",model$beta[s,1],"\n",
      "True beta 2:",model$beta[s,2],"\n",
      "Fitted:\n",
      "Beta 1:",summ$coeff[2,1:2],"\n",
      "Beta 2:",summ$coeff[3,1:2],"\n",
      "Dispersion",summ$dispersion,"\n",
      "Degrees of freedom",summ$df[2],"\n",
      "Correlation",summ$correlation[3,2],"\n\n",     
      file=outfile,append=TRUE)
}
return(NULL) # exists only for side-effect of plotting
}