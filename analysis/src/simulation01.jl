include("./common.jl")

DATAFILE = joinpath(DATADIR, "2015_01_28.txt")

data = readtable(DATAFILE, separator = '\t')
