include("./data.jl")

using CellwiseEditingDifferentiation
using Color
using Gadfly

for i in 1:nsites
  colors = distinguishable_colors(ncells[i])

  layers = Layer[]

  pchain = readdlm(joinpath(OUTDIR, @sprintf("pchain_%s_site%02d.txt", string(simulationid), i)))

  for j in 1:ncells[i]
    push!(layers, layer(
      x=pchain[:, j],
      Stat.histogram(bincount=50, density=true),
      Geom.line,
      Theme(default_color=colors[j])
    )[1])
  end

  vplot = plot(
    layers,
    Guide.xlabel("p<sub>$i</sub>"),
    Guide.title("Histograms of marginal p<sub>$i</sub>"),
    Guide.manual_color_key("Distribution", ["cell "*string(k) for k in 1:length(colors)], [c for c in colors])
  )

  draw(PDF(joinpath(OUTDIR, @sprintf("marginal_phist_%s_site%02d.pdf", string(simulationid), i)), 4inch, 3inch), vplot)
end
