include("./common.jl")

v = readdlm(joinpath(OUTDIR, "simulation01_v.txt"))

vplot = plot(x=v,
             Geom.histogram(bincount=50),
             Guide.xlabel("v"),
             Guide.title("Variance histogram"))

draw(PDF(joinpath(OUTDIR, "simulation01_v_hist.pdf"), 4inch, 3inch), vplot)
