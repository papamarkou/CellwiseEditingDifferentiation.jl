include("./data.jl")

using CellwiseEditingDifferentiation
# using Color
using Gadfly

hyperpars = Dict{Symbol, Any}(:λ=>fill(300.0, nsites))

prior = Dict{Symbol, Any}(:v=>Function[(m::Float64, v::Float64, a::Float64, b::Float64)->
  vpcprior(m, v, a, b, hyperpars[:λ][i]) for i in 1:nsites])

target = Dict{Symbol, Any}(:v=>Array(Function, nsites))
for i in 1:nsites
  target[:v][i] = function (v::Float64)
    a, b = beta_pars_from_mv(data[:m][i], v)
    exp(prior[:v][i](data[:m][i], v, a, b))
  end
end

support = Dict{Symbol, Any}(:v=>Any[0.0001:0.0001:0.2458, 0.0001:0.0001:0.2476])

vprior_xmin = fill(0., nsites)
vprior_xmax = fill(0.03, nsites)
vprior_ymin = fill(0., nsites)
vprior_ymax = fill(200., nsites)

for i in 1:nsites
  # colors = distinguishable_colors(2)
  colors = Dict{Symbol, Color.RGB{Float64}}(:posterior=>color("red"), :prior=>color("blue"))

  layers = Layer[]

  vpdf = Float64[target[:v][i](x) for x in support[:v][i]]
  c = quadgk(target[:v][i], support[:v][i][1], support[:v][i][end])

  push!(layers, layer(
    x=collect(support[:v][i]),
    y=vpdf/c[1],
    Geom.line,
    Theme(default_color=colors[:prior])
  )[1])

  push!(layers, layer(
    x=vec(readdlm(joinpath(OUTDIR, @sprintf("vchain_%s_site%02d.txt", string(simulationid), i)))),
    Stat.histogram(bincount=50, density=true),
    Geom.line,
    Theme(default_color=colors[:posterior])
  )[1])

  vplot = plot(
    layers,
    Guide.xlabel("v<sub>$i</sub>"),
    Guide.title("Histogram of v<sub>$i</sub>"),
    Guide.manual_color_key("Distribution", [string(k) for k in keys(colors)], [c for c in values(colors)]),
    Coord.Cartesian(xmin=vprior_xmin[i], xmax=vprior_xmax[i], ymin=vprior_ymin[i], ymax=vprior_ymax[i])
  )

  draw(PDF(joinpath(OUTDIR, @sprintf("vhist_%s_site%02d.pdf", string(simulationid), i)), 4inch, 3inch), vplot)
end
