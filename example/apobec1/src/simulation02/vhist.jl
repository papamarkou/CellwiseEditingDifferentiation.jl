include("./data.jl")

# using Color
using Gadfly

for i in 1:nsites
  # colors = distinguishable_colors(2)
  colors = Dict{Symbol, Color.RGB{Float64}}(:posterior=>color("red"), :prior=>color("blue"))

  layers = Layer[]

  push!(layers, layer(
    x=vec(readdlm(joinpath(OUTDIR, @sprintf("vchain_%s_site%02d.txt", string(simulationid), i)))),
    Stat.histogram(bincount=50, density=true),
    Geom.line,
    Theme(default_color=colors[:posterior])
  )[1])

  push!(layers, layer(
    x=vec(readdlm(joinpath(OUTDIR, "chain_simulation02_site01.txt"))),
    Stat.histogram(bincount=50, density=true),
    Geom.line,
    Theme(default_color=colors[:prior])
  )[1])

  vplot = plot(
    layers,
    Guide.xlabel("v<sub>$i</sub>"),
    Guide.title("Histogram of v<sub>$i</sub>"),
    Guide.manual_color_key("Distribution", [string(k) for k in keys(colors)], [c for c in values(colors)])
  )

  draw(PDF(joinpath(OUTDIR, @sprintf("vhist_%s_site%02d.pdf", string(simulationid), i)), 4inch, 3inch), vplot)
end
