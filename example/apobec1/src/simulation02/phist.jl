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

include("./data.jl")

# using Color
using Gadfly

for i in 1:nsites
  # colors = distinguishable_colors(2)
  colors = Dict{Symbol, Color.RGB{Float64}}(:posterior=>color("red"), :prior=>color("blue"))

  np = 10000

  layers = Layer[]

  vposterior = vec(readdlm(joinpath(OUTDIR, @sprintf("vchain_%s_site%02d.txt", string(simulationid), i))))

  pposterior = Array(Float64, np)

  for i in 1:np
    a, b = beta_pars_from_mv(m, rand(v))
    pposterior[i] = rand(Beta(a, b))
  end

  println("Mean of pposterior = $(mean(pposterior))")
  println("Var of pposterior = $(var(pposterior))")

  push!(layers, layer(
    x=pposterior,
    Stat.histogram(bincount=50, density=true),
    Geom.line,
    Theme(default_color=colors[:posterior])
  )[1])

  vprior = vec(readdlm(joinpath(OUTDIR, @sprintf("vchain_%s_site%02d.txt", string(simulationid), i))))

  pprior = Array(Float64, np)

  for i in 1:np
    a, b = beta_pars_from_mv(m, rand(v))
    pprior[i] = rand(Beta(a, b))
  end

  println("Mean of pprior = $(mean(pprior))")
  println("Var of pprior = $(var(pprior))")

  push!(layers, layer(
    x=pprior,
    Stat.histogram(bincount=50, density=true),
    Geom.line,
    Theme(default_color=colors[:prior])
  )[1])

  vplot = plot(
    layers,
    Guide.xlabel("p<sub>$i</sub>"),
    Guide.title("Histogram of p<sub>$i</sub>"),
    Guide.manual_color_key("Distribution", [string(k) for k in keys(colors)], [c for c in values(colors)])
  )

  draw(PDF(joinpath(OUTDIR, @sprintf("phist_%s_site%02d.pdf", string(simulationid), i)), 4inch, 3inch), vplot)
end
