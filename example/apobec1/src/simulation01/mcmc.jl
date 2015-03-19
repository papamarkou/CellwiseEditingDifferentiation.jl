include("./artificial_data.jl")

using Distributions
using CellwiseEditingDifferentiation
using Lora

### Targets

hyperpars = Dict{Symbol, Any}(:λ=>fill(10.0, nsites))

prior = Dict{Symbol, Any}(:v=>Function[(m::Float64, v::Float64, a::Float64, b::Float64)->
                          vpcprior(m, v, a, b, hyperpars[:λ][i]) for i in 1:nsites])

target = Dict{Symbol, Any}(:w=>Array(Function, nsites), :p=>Dict{Int, Vector{Function}}())

for i in 1:nsites
  target[:w][i] = function (w::Vector{Float64}, p::Vector{Float64})
    v = data[:m][i]*(1-data[:m][i])*inv_logit(w[1])
    a, b = beta_pars_from_mv(data[:m][i], v)
    sum([logpdf(Beta(a, b), p[m]) for m in 1:ncells[i]])+prior[:v][i](data[:m][i], v, a, b)+logdvdw(w[1], data[:m][i])
  end

  target[:p][i] = Array(Function, ncells[i])
  for j in 1:ncells[i]
    target[:p][i][j] = function (w::Float64)
      a, b = beta_pars_from_mv(data[:m][i], data[:m][i]*(1-data[:m][i])*inv_logit(w))
      Beta(a+data[:edited][i, cells[i][j]], b+data[:coverage][i, cells[i][j]]-data[:edited][i, cells[i][j]])
    end
  end
end

### Initial conditions

init = Dict{Symbol, Any}(:w=>Float64[0. for i in 1:nsites])
init[:p] = Dict{Int, Vector{Float64}}()
for i in 1:nsites
  init[:p][i] = Array(Float64, ncells[i])
  for j in 1:ncells[i]
    init[:p][i][j] = rand(target[:p][i][j](init[:w][i]))
  end
end

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

runner = Dict{Symbol, Any}(:w=>[SerialMC(burnin=0, nsteps=1, thinning=1) for i in 1:nsites])

### Jobs

job = Dict{Symbol, Any}(:w=>Array(Function, nsites), :p=>Dict{Int, Vector{Function}}())

for i in 1:nsites
  job[:w][i] = function (p::Vector{Float64}, init::Vector{Float64})
    MCJob(modeller[:w][i](p, init), sampler[:w][i], runner[:w][i])
  end

  job[:p][i] = Array(Function, ncells[i])
  for j in 1:ncells[i]
    job[:p][i][j] = function (w::Float64)
      DistJob(target[:p][i][j](w))
    end
  end
end

job[:gibbs] = SerialMC(burnin=1000, nsteps=11000, thinning=1)

### Run Gibbs sampler

mcchain = gibbs(data, init, job)

### Write output

for i in 1:nsites
  VOUTFILE = joinpath(OUTDIR, @sprintf("vchain_%s_site%02d.txt", string(simulationid), i))
  writedlm(VOUTFILE,  map(w->data[:m][i]*(1-data[:m][i])*inv_logit(w), mcchain[:w][:, i]), ' ')

  POUTFILE = joinpath(OUTDIR, @sprintf("pchain_%s_site%02d.txt", string(simulationid), i))
  writedlm(POUTFILE, mcchain[:p][i], ' ')
end
