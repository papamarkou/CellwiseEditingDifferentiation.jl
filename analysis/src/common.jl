using DataFrames
using Distributions
using Lora

DATADIR = "../../data"
OUTDIR = "../output"

include("./logit.jl")
include("./beta_moments.jl")
include("./gibbs.jl")