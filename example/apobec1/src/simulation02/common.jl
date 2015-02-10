include("../common.jl")

DATADIR = "../../data"
OUTROOT = "../output"

simulationid = :simulation02

OUTDIR = joinpath(OUTROOT, string(simulationid))

### Load data

DATAFILE = joinpath(DATADIR, "2015_01_28.txt")

df = readtable(DATAFILE, separator = '\t')
df = df[df[:coord] .== 122152902, :]

data = Dict{Symbol, Any}(:m=>df[1, :Bulk_er],
                         :coverage=>float(array(df[1, 5:2:51])),
                         :rate=>float(array(df[1, 6:2:52])))
data[:edited] = float(data[:rate].*data[:coverage])
nsites, ncells = size(data[:coverage])
