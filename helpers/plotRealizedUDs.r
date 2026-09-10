# Generate figures

###

conditional_ud <- UD_out$conditional_ud
joint_ud <- UD_out$joint_ud
marginal_ud <- UD_out$marginal_ud
hh <- UD_out$hh

header_naive <- "(a) Independent transition rates"
header_novel <- "(b) Joint transition rates"

#
# KDE plot - stacked
#
if (PDFs) pdf(paste0(destination,scenario,"_KDEs_beta",paste0(betaCoeff,collapse="_"),"_sigma",paste0(sig,collapse="_"),".pdf"),height=10,width=7)
par(mar=c(4,5,4,3))
layout(as.matrix(1:2))#,respect=TRUE)

# Naive KDE
initial_plots(model,controls,UD_out)
out_pooled_naive <- display_pooled(model$S,controls,appear,
                                   outputs_naive,header_naive)

# Novel KDE
initial_plots(model,controls,UD_out)
out_pooled_novel <- display_pooled(model$S,controls,appear,
                                   outputs_novel,header_novel)

if (PDFs) dev.off()

#
# Cell probabilities
#
if (PDFs) pdf(paste0(destination,scenario,"_Bins_beta",paste0(betaCoeff,collapse="_"),"_sigma",paste0(sig,collapse="_"),".pdf"),height=9,width=6)
par(mar=c(4,5,4,3))
layout(matrix(1:6,byrow=FALSE,ncol=2),respect=TRUE)

# Naive bins
bin_out_naive <- bins(model$S,controls,appear,UD_out,out_pooled_naive,header_naive)

if(GLMs & model$K==2)
{
  outfile <- paste0(destination,scenario,"GLMs_alt.txt")
  glms_naive <- fit(model,habitat,bin_out_naive)
  outfit(model,glms_naive,outfile,scenario,header_naive)
}

# Novel bins
bin_out_novel <- bins(model$S,controls,appear,UD_out,out_pooled_novel,header_novel)


if(GLMs & model$K==2)
{
  glms_novel <- fit(model,habitat,bin_out_novel)
  outfit(model,glms_novel,outfile,scenario,header_novel)
}
if (PDFs) dev.off()

