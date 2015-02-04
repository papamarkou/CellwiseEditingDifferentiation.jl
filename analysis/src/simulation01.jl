include("./common.jl")

DATAFILE = joinpath(DATADIR, "2015_01_28.txt")

df = readtable(DATAFILE, separator = '\t')
df = df[df[:coord] .== 122152902, :]

data = Dict{Symbol, Any}(:m=>df[1, :Bulk_er],
                         :coverage=>float(array(df[1, 5:2:51])),
                         :rate=>float(array(df[1, 6:2:52])))
data[:edited] = float(data[:rate].*data[:coverage])
nsites, ncells = size(data[:coverage])

gibbs_prior = Dict{Symbol, Any}(:w=>Distribution[Normal() for i in 1:nsites])

gibbs_init = Dict{Symbol, Any}(:v=>Float64[rand(gibbs_prior[:v][i]) for i in 1:nsites])
p = Array(Float64, nsites, ncells)
for i in 1:nsites
  a, b = beta_pars_from_mv(data[:m], gibbs_init[:v][i])
  for j in 1:ncells
    p[i, :] = rand(Beta(a+data[:edited][i, j], b+data[:coverage][i, j]-data[:edited][i, j]))
  end
end
gibbs_init[:p] = p

gibbs_runner = Dict{Symbol, Any}(:burnin=>100,
                                 :nsteps=>1000,
                                 :thinning=>1)

inner_runner = Dict{Symbol, Any}(:burnin=>100,
                                 :nsteps=>1000,
                                 :thinning=>1)

mcchain = gibbs(data, gibbs_prior, gibbs_init, gibbs_runner, metropolis_prior, inner_runner)

VOUTFILE = joinpath(OUTDIR, "simulation01_v.txt")
writedlm(VOUTFILE, mcchain[:v], ' ')

POUTFILES = [joinpath(OUTDIR, @sprintf("simulation01d_p_cell%02d.txt", j)) for j in 1:ncells]
for j in 1:ncells
  writedlm(POUTFILES[j], mcchain[:p][:, :, j], ' ')
end
