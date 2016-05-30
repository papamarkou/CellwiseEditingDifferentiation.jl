simulationid = :simulation01

OUTDIR = joinpath("../../output", string(simulationid))

nsites = 1

data = Dict{Symbol, Any}(:m=>0.5,
                         :coverage=>fill(10.0, 20)',
                         :edited=>rand([1., 2., 3., 8., 9.], 20)')
data[:rate] = float(data[:edited]./data[:coverage])

cells = Dict{Int, Vector{Int}}()
for i in 1:nsites
  cells[i] = find(data[:coverage][i, :] .!= 0)
end

ncells = Dict{Int, Int}()
for i in 1:nsites
  ncells[i] = length(cells[i])
end
