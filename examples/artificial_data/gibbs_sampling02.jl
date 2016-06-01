using CellwiseEditingDifferentiation
using Distributions
using Lora

ncells = 20

λ = Hyperparameter(:λ)

μ = Data(:μ)
coverage = Data(:coverage)
edited = Data(:edited)
rate = Data(:rate)

function wlogtarget(w::Float64, s::Vector)
  v = s[2]*(1-s[2])*inv_logit(w)
  a, b = beta_pars_from_mv(s[2], v)
  mapreduce(x -> logpdf(Beta(a, b), x), +, s[7:(6+ncells)])+vpcprior(s[2], v, a, b, s[1])+logdvdw(w, s[2])
end

w = BasicContUnvParameter(:w, logtarget=wlogtarget, nkeys=6+ncells)

psetpdf = Array(Function, ncells)
for i in 1:ncells
  psetpdf[i] = function (p::Float64, s::Vector)
    a, b = beta_pars_from_mv(s[2], s[2]*(1-s[2])*inv_logit(s[6]))
    Beta(a+s[4][i], b+s[3][i]-s[4][i])
  end
end

p = Array(BasicContUnvParameter, ncells)
for i in 1:ncells
  p[i] = BasicContUnvParameter(Symbol('p', i), setpdf=psetpdf[i], nkeys=6+ncells)
end

model = GenericModel([λ, μ, coverage, edited, rate, w, p...], isindexed=false)

init = Dict(
  :λ => 10.,
  :μ => 0.5,
  :coverage => fill(10.0, 20),
  :edited => rand([1., 2., 3., 8., 9.], 20),
  :w => 0.
)
init[:rate] = float(init[:edited]./init[:coverage])
for i in 1:ncells
  init[Symbol('p', i)] =
    rand(psetpdf[i](NaN, [[init[k] for k in [:λ, :μ, :coverage, :edited, :rate, :w]]; fill(NaN, ncells)]))
end
  
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

outopts = fill(Dict(:monitor => [:value]), ncells+1)

gibbsjob = GibbsJob(model, Dict(:w => wjob), mcrange, init)
  
@time run(gibbsjob)

output(gibbsjob)

chains = Dict(gibbsjob)

v = map(x -> init[:μ]*(1-init[:μ])*inv_logit(x), chains[:w].value)
  
for i in 1:ncells
  writecsv(string("pchain", i, ".csv"), chains[Symbol('p', i)].value)  
end
writecsv("vchain.csv", v)
