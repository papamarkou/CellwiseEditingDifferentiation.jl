using DataFrames

simulationid = :simulation01

DATADIR = "../../data"
OUTDIR = joinpath("../../output", string(simulationid))

### Load data

DATAFILE = joinpath(DATADIR, "apobec1.txt")

df = readtable(DATAFILE, separator = '\t')

sites = DataFrame(chr = ["chr10", "chr11"], coord = [57515985, 100476998])
# sites = DataFrame(chr = ["chr2", "chr16"], coord = [122152902, 84954513])
# sites = DataFrame(chr = ["chrX", "chr4", "chr5"], coord = [109162174, 6395365, 104440840])
nsites = nrow(sites)

df = df[intersect(findin(df[:chr], sites[:chr]), findin(df[:coord], sites[:coord])), :]

data = Dict{Symbol, Any}(:m=>df[:, :Bulk_er],
                         :coverage=>float(array(df[:, 5:2:51])),
                         :rate=>float(array(df[:, 6:2:52])))
data[:edited] = float(data[:rate].*data[:coverage])

cells = Dict{Int, Vector{Int}}()
for i in 1:nsites
  cells[i] = find(data[:coverage][i, :] .!= 0)
end

ncells = Dict{Int, Int}()
for i in 1:nsites
  ncells[i] = length(cells[i])
end
