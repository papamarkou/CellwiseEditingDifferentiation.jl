using DataFrames
using Distributions
using Gadfly
using Lora

DATADIR = "../../data"
OUTROOT = "../output"

include("./logit.jl")
include("./beta_moments.jl")
include("./DistJob.jl")
include("./gibbs.jl")
include("./vprior.jl")
