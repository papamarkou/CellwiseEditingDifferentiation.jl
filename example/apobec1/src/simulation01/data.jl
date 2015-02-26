using DataFrames

simulationid = :simulation01

DATADIR = "../../data"
OUTDIR = joinpath("../../output", string(simulationid))

### Load data

DATAFILE = joinpath(DATADIR, "apobec1.txt")

sites = DataFrame(chr = ["chr2", "chr16"], coord = [122152902, 84954513])
# sites = DataFrame(chr = ["chrX", "chr4", "chr5"], coord = [109162174, 6395365, 104440840])
nsites = nrow(sites)

df = readtable(DATAFILE, separator = '\t')
df = df[intersect(findin(df[:chr], sites[:chr]), findin(df[:coord], sites[:coord])), :]

data = Dict{Symbol, Any}(:m=>df[:, :Bulk_er],
                         :coverage=>float(array(df[:, 5:2:51])),
                         :rate=>float(array(df[:, 6:2:52])))
data[:edited] = float(data[:rate].*data[:coverage])
ncells = size(data[:coverage], 2)
