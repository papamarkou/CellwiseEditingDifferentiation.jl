include("./common.jl")

v = vec(readdlm(joinpath(OUTDIR, "simulation02_v.txt")))

vplot = plot(x=v,
             Geom.histogram(bincount=50),
             Guide.xlabel("v"),
             Guide.title("Variance histogram"))

draw(PDF(joinpath(OUTDIR, "simulation02_vhist.pdf"), 4inch, 3inch), vplot)
