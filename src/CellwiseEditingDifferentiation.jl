module CellwiseEditingDifferentiation

using Distributions
using Lora

import Base:
  run

export
  logit,
  inv_logit,
  beta_a_from_mv,
  beta_b_from_mv,
  beta_pars_from_mv,
  vprior,
  logdvdw,
  wprior,
  DistJob,
  Job,
  run,
  gibbs

include("logit.jl")
include("beta_moments.jl")
include("vprior.jl")
include("DistJob.jl")
include("gibbs.jl")

end
