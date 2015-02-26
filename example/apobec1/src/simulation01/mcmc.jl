include("./data.jl")

using Distributions
using CellwiseEditingDifferentiation
using Lora

### Targets

hyperpars = Dict{Symbol, Any}(:λ=>fill(300.0, nsites))

prior = Dict{Symbol, Any}(:v=>Function[(m::Float64, v::Float64, a::Float64, b::Float64)->
                          vpcprior(m, v, a, b, hyperpars[:λ][i]) for i in 1:nsites])

target = Dict{Symbol, Any}(:w=>Array(Function, nsites), :p=>Array(Function, nsites, ncells))

for i in 1:nsites
  target[:w][i] = function (w::Vector{Float64}, p::Vector{Float64})
    v = data[:m][i]*(1-data[:m][i])*inv_logit(w[1])
    a, b = beta_pars_from_mv(data[:m][i], v)
    sum([logpdf(Beta(a+data[:edited][i, m], b+data[:coverage][i, m]-data[:edited][i, m]), p[m]) for m in 1:ncells])+
      2*logdvdw(w[1], data[:m][i])+
      prior[:v][i](data[:m][i], v, a, b)
  end

  for j in 1:ncells
    target[:p][i, j] = function (w::Float64)
      a, b = beta_pars_from_mv(data[:m][i], data[:m][i]*(1-data[:m][i])*inv_logit(w))
      Beta(a+data[:edited][i, j], b+data[:coverage][i, j]-data[:edited][i, j])
    end
  end
end

### Initial conditions

init = Dict{Symbol, Any}(:w=>Float64[0. for i in 1:nsites])
init[:p] = Float64[rand(target[:p][i, j](init[:w][i])) for i in 1:nsites, j in 1:ncells]

### Modellers

modeller = Dict{Symbol, Any}(:w=>Array(Function, nsites))
for i in 1:nsites
  modeller[:w][i] = function (p::Vector{Float64}, init::Vector{Float64})
    model(w->target[:w][i](w, p), init=init)
  end
end

### Samplers

sampler = Dict{Symbol, Any}(:w=>[RAM() for i in 1:nsites])

### Runners

runner = Dict{Symbol, Any}(:w=>[SerialMC(burnin=1000, nsteps=11000, thinning=1) for i in 1:nsites])

### Jobs

job = Dict{Symbol, Any}(:w=>Array(Function, nsites), :p=>Array(Function, nsites, ncells))

for i in 1:nsites
  job[:w][i] = function (p::Vector{Float64}, init::Vector{Float64})
    MCJob(modeller[:w][i](p, init), sampler[:w][i], runner[:w][i])
  end

  for j in 1:ncells
    job[:p][i, j] = function (w::Float64)
      DistJob(target[:p][i, j](w))
    end
  end
end

job[:gibbs] = SerialMC(burnin=1000, nsteps=11000, thinning=1)

### Run Gibbs sampler

mcchain = gibbs(data, init, job)

### Write output

for i in 1:nsites
  VOUTFILE = joinpath(OUTDIR, @sprintf("vchain_%s_site%02d.txt", string(simulationid), i))
  writedlm(VOUTFILE,  map(w->0.25*inv_logit(w), mcchain[:w][:, i]), ' ')

  for j in 1:ncells
    POUTFILE = joinpath(OUTDIR, @sprintf("pchain_%s_site%02d_cell%02d.txt", string(simulationid), i, j))
    writedlm(POUTFILE, mcchain[:p][:, i, j], ' ')
  end
end
