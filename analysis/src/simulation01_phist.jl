include("./common.jl")

### Read Beta mean from data

DATAFILE = joinpath(DATADIR, "2015_01_28.txt")

df = readtable(DATAFILE, separator = '\t')
df = df[df[:coord] .== 122152902, :]

m = df[1, :Bulk_er]

### Read Beta variance from MCMC simulation

v = vec(readdlm(joinpath(OUTDIR, "simulation01_v.txt")))

### Draw samples from Beta distribution of editing rate

np = 10000

p = Array(Float64, np)

for i in 1:np
  a, b = beta_pars_from_mv(m, rand(v))
  p[i] = rand(Beta(a, b))
end

println("Mean of p = $(mean(p))")
println("Var of p = $(var(p))")

vplot = plot(x=p,
             Geom.histogram(bincount=50),
             Guide.xlabel("v"),
             Guide.title("Editing rate histogram"))

draw(PDF(joinpath(OUTDIR, "simulation01_phist.pdf"), 4inch, 3inch), vplot)
