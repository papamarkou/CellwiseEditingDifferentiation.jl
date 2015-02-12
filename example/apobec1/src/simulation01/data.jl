using DataFrames

simulationid = :simulation02

DATADIR = "../../data"
OUTDIR = joinpath("../../output", string(simulationid))

### Load data

DATAFILE = joinpath(DATADIR, "apobec1.txt")

sites = [122152902, 84954513]
nsites = length(sites)

df = readtable(DATAFILE, separator = '\t')
df = df[findin(df[:coord], sites), :]

data = Dict{Symbol, Any}(:m=>df[:, :Bulk_er],
                         :coverage=>float(array(df[:, 5:2:51])),
                         :rate=>float(array(df[:, 6:2:52])))
data[:edited] = float(data[:rate].*data[:coverage])
ncells = size(data[:coverage], 2)
