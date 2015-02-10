include("./common.jl")

for i in 1:nsites
  VOUTFILE = joinpath(OUTDIR, @sprintf("vchain_%s_site%02d.txt", string(simulationid), i))

  vplot = plot(x=vec(readdlm(VOUTFILE)),
               Geom.histogram(bincount=50),
               Guide.xlabel("v"),
               Guide.title("Variance histogram"))

  draw(PDF(joinpath(OUTDIR, "simulation02_vhist.pdf"), 4inch, 3inch), vplot)
end
