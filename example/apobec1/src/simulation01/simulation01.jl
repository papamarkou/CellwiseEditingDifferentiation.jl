using CellwiseEditingDifferentiation
using Distributions
using Lora

import Base: rand

immutable MvBeta <: ContinuousMultivariateDistribution
  betas::Vector{Beta}
end

Base.rand(d::MvBeta) = [rand(d.betas[i]) for i in 1:length(d.betas)]

ncells = 20

λ = Hyperparameter(:λ)

μ = Data(:μ)
coverage = Data(:coverage)
edited = Data(:edited)
rate = Data(:rate)

w = BasicContUnvParameter(
  :w,
  logtarget = function (w::Float64, v::Vector)
    θ = v[2]*(1-v[2])*inv_logit(w)
    a, b = beta_pars_from_mv(v[2], θ)
    mapreduce(x -> logpdf(Beta(a, b), x), +, v[7])+vpcprior(v[2], θ, a, b, v[1])+logdvdw(w, v[2])
  end,
  nkeys=7
)

p = BasicContMuvParameter(
  :p,
  setpdf = function (p::Vector{Float64}, v::Vector)
    a, b = beta_pars_from_mv(v[2], v[2]*(1-v[2])*inv_logit(v[6]))
    MvBeta([Beta(a+v[4][i], b+v[3][i]-v[4][i]) for i in 1:ncells])
  end,
  nkeys=7
)

model = GenericModel(
  [λ, μ, coverage, edited, rate, w, p],
  [λ w; μ w; μ p; coverage p; edited p; w p],
  isindexed=false
)

v0 = Dict(
  :λ => 10.,
  :μ => 0.5,
  :coverage => fill(10.0, 20),
  :edited => rand([1., 2., 3., 8., 9.], 20),
  :w => 0.,
  :p => fill(0.1, ncells)
)
v0[:rate] = float(v0[:edited]./v0[:coverage])
  
wjob = BasicMCJob(
  model,
  RAM(),
  BasicMCRange(nsteps=100, burnin=10),
  v0;
  pindex=6,
  resetpstate=false,
  verbose=false
)

mcrange = BasicMCRange(nsteps=1100, burnin=100)

outopts = [Dict(:monitor => [:value]), Dict(:monitor => [:value])]

gibbsjob = GibbsJob(model, Dict(:w => wjob), mcrange, v0)
  
@time run(gibbsjob)

output(gibbsjob)

chains = Dict(gibbsjob)

θ = map(x -> v0[:μ]*(1-v0[:μ])*inv_logit(x), chains[:w].value)

