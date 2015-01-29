include("./common.jl")

DATAFILE = joinpath(DATADIR, "2015_01_28.txt")

df = readtable(DATAFILE, separator = '\t')
df = df[df[:coord] .== 122152902, :]

data = Dict{Symbol, Any}(:m=>df[1, :Bulk_er],
                         :coverage=>float(array(df[1, 5:2:51])),
                         :rate=>float(array(df[1, 6:2:52])))
data[:edited] = float(data[:rate].*data[:coverage])
nsites, ncells = size(data[:coverage])

gibbs_prior = Dict{Symbol, Any}(:v=>Distribution[Beta(1.0, 1.0) for i in 1:nsites])

gibbs_init = Dict{Symbol, Any}(:v=>Float64[rand(gibbs_prior[:v][i]) for i in 1:nsites])
p = Array(Float64, nsites, ncells)
for i in 1:nsites
  a, b = beta_pars_from_mv(data[:m], gibbs_init[:v][i])
  for j in 1:ncells
    p[i, :] = rand(Beta(a+data[:edited][i, j], b+data[:coverage][i, j]-data[:edited][i, j]))
  end
end
gibbs_init[:p] = p

gibbs_runner = Dict{Symbol, Int}(:burnin=>10,
                                 :nsteps=>100,
                                 :thinning=>1)

metropolis_prior = Dict{Symbol, Distribution}(:Σ=>fill(eye(1), nsites))

metropolis_runner = Dict{Symbol, Int}(:burnin=>1000,
                                      :nsteps=>10000,
                                      :thinning=>1)

v, p = metropolis_within_gibbs(data, gibbs_prior, gibbs_init, gibbs_runner, metropolis_prior, metropolis_runner)
