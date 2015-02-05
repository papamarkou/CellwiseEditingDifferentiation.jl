using DataFrames
using Distributions
using Gadfly
using Lora

DATADIR = "../../data"
OUTDIR = "../output"

include("./logit.jl")
include("./beta_moments.jl")
include("./DistJob.jl")
include("./gibbs.jl")
