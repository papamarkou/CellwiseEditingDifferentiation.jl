module CellwiseEditingDifferentiation

using Distributions
using Lora

export
  logit,
  inv_logit,
  beta_a_from_mv,
  beta_b_from_mv,
  beta_pars_from_mv,
  vpcprior,
  logdvdw,
  wpcprior

include("logit.jl")
include("beta_moments.jl")
include("vpcprior.jl")

end
