include("./data.jl")

using Distributions
using CellwiseEditingDifferentiation
using Lora

### Targets

hyperpars = Dict{Symbol, Any}(:λ=>fill(300.0, nsites))

prior = Dict{Symbol, Any}(:w=>Function[(m::Float64, v::Float64, w::Float64, a::Float64, b::Float64)->
                          wprior(m, v, w, a, b, hyperpars[:λ][i]) for i in 1:nsites])

target = Dict{Symbol, Any}(:w=>Array(Function, nsites))

for i in 1:nsites
  target[:w][i] = function (w::Vector{Float64})
    v = 0.25*inv_logit(w[1])
    a, b = beta_pars_from_mv(data[:m][i], v)
    prior[:w][i](data[:m][i], v, w[1], a, b)
  end
end

### Initial conditions

init = Dict{Symbol, Any}(:w=>Float64[0. for i in 1:nsites])

### Modellers

modeller = Dict{Symbol, Any}(:w=>Array(Function, nsites))
for i in 1:nsites
  modeller[:w][i] = function (init::Vector{Float64})
    model(w->target[:w][i](w), init=init)
  end
end

### Samplers

sampler = Dict{Symbol, Any}(:w=>[RAM() for i in 1:nsites])

### Runners

runner = Dict{Symbol, Any}(:w=>[SerialMC(burnin=10, nsteps=110, thinning=1) for i in 1:nsites])

### Jobs

job = Dict{Symbol, Any}(:w=>Array(Function, nsites))

for i in 1:nsites
  job[:w][i] = function (init::Vector{Float64})
    MCJob(modeller[:w][i](init), sampler[:w][i], runner[:w][i])
  end
end

### Run RAM sampler

for i in 1:nsites
  mcchain = Lora.run(job[:w][i]([init[:w][i]])).samples[:]

  VOUTFILE = joinpath(OUTDIR, @sprintf("vchain_%s_site%02d.txt", string(simulationid), i))

  writedlm(VOUTFILE,  map(w->0.25*inv_logit(w), mcchain), ' ')
end
