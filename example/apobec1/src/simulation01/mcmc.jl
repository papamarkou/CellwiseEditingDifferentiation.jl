include("./common.jl")

### Data

DATAFILE = joinpath(DATADIR, "2015_01_28.txt")

df = readtable(DATAFILE, separator = '\t')
df = df[df[:coord] .== 122152902, :]

data = Dict{Symbol, Any}(:m=>df[1, :Bulk_er],
                         :coverage=>float(array(df[1, 5:2:51])),
                         :rate=>float(array(df[1, 6:2:52])))
data[:edited] = float(data[:rate].*data[:coverage])
nsites, ncells = size(data[:coverage])

### Targets

hyperpars = Dict{Symbol, Any}(:λ=>fill(300.0, nsites))

prior = Dict{Symbol, Any}(:w=>Function[(m::Float64, v::Float64, w::Float64, a::Float64, b::Float64)->
                          wprior(m, v, w, a, b, hyperpars[:λ][i]) for i in 1:nsites])

target = Dict{Symbol, Any}(:w=>Array(Function, nsites), :p=>Array(Function, nsites, ncells))

for i in 1:nsites
  target[:w][i] = function (w::Vector{Float64}, p::Vector{Float64})
    v = 0.25*inv_logit(w[1])
    a, b = beta_pars_from_mv(data[:m], v)
    sum([logpdf(Beta(a+data[:edited][i, m], b+data[:coverage][i, m]-data[:edited][i, m]), p[m]) for m in 1:ncells])+
      logdvdw(w[1])+
      prior[:w][i](data[:m], v, w[1], a, b)
  end

  for j in 1:ncells
    target[:p][i, j] = function (w::Float64)
      a, b = beta_pars_from_mv(data[:m], 0.25*inv_logit(w))
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

VOUTFILE = joinpath(OUTDIR, "simulation01_v.txt")
writedlm(VOUTFILE,  map(w->0.25*inv_logit(w), mcchain[:w]), ' ')

POUTFILES = [joinpath(OUTDIR, @sprintf("simulation01d_p_cell%02d.txt", j)) for j in 1:ncells]
for j in 1:ncells
  writedlm(POUTFILES[j], mcchain[:p][:, :, j], ' ')
end
