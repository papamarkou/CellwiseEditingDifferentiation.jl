include("./data.jl")

using CellwiseEditingDifferentiation
# using Color
using Distributions
using Gadfly
using Lora

hyperpars = Dict{Symbol, Any}(:λ=>fill(300.0, nsites))

prior = Dict{Symbol, Any}(:v=>Function[(m::Float64, v::Float64, a::Float64, b::Float64)->
  vpcprior(m, v, a, b, hyperpars[:λ][i]) for i in 1:nsites])

prior_target = Dict{Symbol, Any}(:w=>Array(Function, nsites))
for i in 1:nsites
  prior_target[:w][i] = function (w::Vector{Float64})
    v = data[:m][i]*(1-data[:m][i])*inv_logit(w[1])
    a, b = beta_pars_from_mv(data[:m][i], v)
    prior[:v][i](data[:m][i], v, a, b)+logdvdw(w[1], data[:m][i])
  end
end

prior_init = Dict{Symbol, Any}(:w=>fill(0.0, nsites))

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

  mcmodel = model(prior_target[:w][i], init=[prior_init[:w][i]])
  mcsampler = RAM()
  mcrunner = SerialMC(nsteps=int(1.01*np), burnin=int(0.1*np))
  mcchain = run(mcmodel, mcsampler, mcrunner)
  mcchain = map(w->data[:m][i]*(1-data[:m][i])*inv_logit(w), vec(mcchain.samples))

  pprior = Array(Float64, np)

  for j in 1:ncells
    a, b = beta_pars_from_mv(data[:m][i], rand(mcchain))
    pprior[j] = rand(Beta(a+data[:edited][i, j], b+data[:coverage][i, j]-data[:edited][i, j]))
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
