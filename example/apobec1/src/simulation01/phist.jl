include("./data.jl")

using CellwiseEditingDifferentiation
# using Color
using Distributions
using Gadfly

hyperpars = Dict{Symbol, Any}(:λ=>fill(300.0, nsites))

prior = Dict{Symbol, Any}(:v=>Function[(m::Float64, v::Float64, a::Float64, b::Float64)->
  vpcprior(m, v, a, b, hyperpars[:λ][i]) for i in 1:nsites])

target = Dict{Symbol, Any}(:v=>Array(Function, nsites))
for i in 1:nsites
  target[:v][i] = function (v::Float64)
    a, b = beta_pars_from_mv(data[:m][i], v)
    try
      exp(prior[:v][i](data[:m][i], v, a, b))
    catch
      0.
    end
  end
end

support = Dict{Symbol, Any}(:v=>Any[0.0001:0.0001:0.2458, 0.0001:0.0001:0.2476])

for i in 1:nsites
  np = 100000

  vposterior = vec(readdlm(joinpath(OUTDIR, @sprintf("vchain_%s_site%02d.txt", string(simulationid), i))))

  pposterior = Array(Float64, np)

  for j in 1:np
    a, b = beta_pars_from_mv(data[:m][i], rand(vposterior))
    pposterior[j] = rand(Beta(a, b))
  end

  println("Mean of pposterior for site $i = $(mean(pposterior))")
  println("Var of pposterior for site $i = $(var(pposterior))")

  u = 0.25*rand(Uniform(0, 1), np)
  weight = map(target[:v][i], u)
  weight = weight/sum(weight)
  index = rand(Categorical(weight), np)
  vprior = u[index]

  pprior = Array(Float64, np)

  for j in 1:np
    a, b = beta_pars_from_mv(data[:m][i], vprior[i])
    pprior[j] = rand(Beta(a, b))
  end

  println("Mean of pprior for site $i = $(mean(pprior))")
  println("Var of pprior for site $i = $(var(pprior))")

  # colors = distinguishable_colors(2)
  colors = Dict{Symbol, Color.RGB{Float64}}(:posterior=>color("red"), :prior=>color("blue"))

  layers = Layer[]

  push!(layers, layer(
    x=pposterior,
    Stat.histogram(bincount=50, density=true),
    Geom.line,
    Theme(default_color=colors[:posterior])
  )[1])

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
