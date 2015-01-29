using Distributions
using DataFrames

DATADIR = "../../data"
OUTDIR = "../output"

include("beta_moments.jl")
include("./metropolis.jl")
include("./metropolis_within_gibbs.jl")
