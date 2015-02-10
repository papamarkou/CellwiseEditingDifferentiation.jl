include("./common.jl")

### Read Beta variance from MCMC simulation

v = vec(readdlm(joinpath(OUTDIR, "simulation02_v.txt")))

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
             Guide.xlabel("p"),
             Guide.title("Editing rate histogram"))

draw(PDF(joinpath(OUTDIR, "simulation02_phist.pdf"), 4inch, 3inch), vplot)
