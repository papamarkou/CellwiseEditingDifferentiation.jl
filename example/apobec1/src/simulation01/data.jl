using DataFrames

simulationid = :simulation01

DATADIR = "../../data"
OUTDIR = joinpath("../../output", string(simulationid))

### Load data

DATAFILE = joinpath(DATADIR, "apobec1.txt")

sites = Dict{ASCIIString, Int}("chr2"=>122152902, "chr16"=>84954513)
nsites = length(sites)

df = readtable(DATAFILE, separator = '\t')
df = df[intersect(findin(df[:chr], keys(sites)), findin(df[:coord], values(sites))), :]

data = Dict{Symbol, Any}(:m=>df[:, :Bulk_er],
                         :coverage=>float(array(df[:, 5:2:51])),
                         :rate=>float(array(df[:, 6:2:52])))
data[:edited] = float(data[:rate].*data[:coverage])
ncells = size(data[:coverage], 2)
