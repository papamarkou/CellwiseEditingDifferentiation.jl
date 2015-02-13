include("./data.jl")

using CellwiseEditingDifferentiation
# using Color
using Distributions
using Gadfly

for i in 1:nsites
  # colors = distinguishable_colors(2)
  colors = Dict{Symbol, Color.RGB{Float64}}(:posterior=>color("red"))

  np = 10000

  layers = Layer[]

  vposterior = vec(readdlm(joinpath(OUTDIR, @sprintf("vchain_%s_site%02d.txt", string(simulationid), i))))

  pposterior = Array(Float64, np)

  for j in 1:np
    a, b = beta_pars_from_mv(data[:m][i], rand(vposterior))
    pposterior[j] = rand(Beta(a, b))
  end

  println("Mean of pposterior = $(mean(pposterior))")
  println("Var of pposterior = $(var(pposterior))")

  push!(layers, layer(
    x=pposterior,
    Stat.histogram(bincount=50, density=true),
    Geom.line,
    Theme(default_color=colors[:posterior])
  )[1])

  vplot = plot(
    layers,
    Guide.xlabel("p<sub>$i</sub>"),
    Guide.title("Histogram of p<sub>$i</sub>"),
    #Guide.manual_color_key("Distribution", [string(k) for k in keys(colors)], [c for c in values(colors)])
  )

  draw(PDF(joinpath(OUTDIR, @sprintf("phist_%s_site%02d.pdf", string(simulationid), i)), 4inch, 3inch), vplot)
end
