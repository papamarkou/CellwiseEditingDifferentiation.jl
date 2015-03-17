simulationid = :simulation01

OUTDIR = joinpath("../../output", string(simulationid))

nsites = 1

# data = Dict{Symbol, Any}(:m=>0.5,
#                          :coverage=>Float64[10 20 30 40],
#                          :edited=>Float64[5 10 15 20])
data = Dict{Symbol, Any}(:m=>0.5,
                         :coverage=>Float64[10 20 30 40],
                         :edited=>Float64[1 2 28 29])
data[:rate] = float(data[:edited]./data[:coverage])

cells = Dict{Int, Vector{Int}}()
for i in 1:nsites
  cells[i] = find(data[:coverage][i, :] .!= 0)
end

ncells = Dict{Int, Int}()
for i in 1:nsites
  ncells[i] = length(cells[i])
end
