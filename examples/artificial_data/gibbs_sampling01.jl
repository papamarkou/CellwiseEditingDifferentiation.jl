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

function wlogtarget(w::Float64, s::Vector)
  v = s[2]*(1-s[2])*inv_logit(w)
  a, b = beta_pars_from_mv(s[2], v)
  mapreduce(x -> logpdf(Beta(a, b), x), +, s[7])+vpcprior(s[2], v, a, b, s[1])+logdvdw(w, s[2])
end

w = BasicContUnvParameter(:w, logtarget=wlogtarget, nkeys=7)

function psetpdf(p::Vector{Float64}, s::Vector)
  a, b = beta_pars_from_mv(s[2], s[2]*(1-s[2])*inv_logit(s[6]))
  MvBeta([Beta(a+s[4][i], b+s[3][i]-s[4][i]) for i in 1:ncells])
end

p = BasicContMuvParameter(:p, setpdf=psetpdf, nkeys=7)

model = GenericModel(
  [λ, μ, coverage, edited, rate, w, p],
  [λ w; μ w; μ p; coverage p; edited p; w p],
  isindexed=false
)

init = Dict(
  :λ => 10.,
  :μ => 0.5,
  :coverage => fill(10.0, 20),
  :edited => rand([1., 2., 3., 8., 9.], 20),
  :w => 0.
)
init[:rate] = float(init[:edited]./init[:coverage])
init[:p] = rand(psetpdf(fill(NaN, ncells), [init[k] for k in [:λ, :μ, :coverage, :edited, :rate, :w]]))
  
wjob = BasicMCJob(
  model,
  RAM(),
  BasicMCRange(nsteps=100, burnin=10),
  init;
  pindex=6,
  resetpstate=false,
  verbose=false
)

mcrange = BasicMCRange(nsteps=1100, burnin=100)

outopts = [Dict(:monitor => [:value]), Dict(:monitor => [:value])]

gibbsjob = GibbsJob(model, Dict(:w => wjob), mcrange, init)
  
@time run(gibbsjob)

output(gibbsjob)

chains = Dict(gibbsjob)

v = map(x -> init[:μ]*(1-init[:μ])*inv_logit(x), chains[:w].value)

writecsv("pchains.csv", chains[:p].value)  
writecsv("vchain.csv", v)
