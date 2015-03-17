### This is not a generically coded Metropolis-within-Gibbs sampler, it is tailored to the problem at hand

function gibbs(data::Dict{Symbol, Any}, init::Dict{Symbol, Any}, job)
  nsites = size(data[:coverage], 1)

  cells = Dict{Int, Vector{Int}}()
  for i in 1:nsites
    cells[i] = find(data[:coverage][i, :] .!= 0)
  end

  ncells = Dict{Int, Int}()
  for i in 1:nsites
    ncells[i] = length(cells[i])
  end

  nmc = length(job[:gibbs].r)

  mcsample = Dict{Symbol, Any}(:w=>init[:w], :p=>init[:p])

  mcchain = Dict{Symbol, Any}(:w=>Array(Float64, nmc, nsites))
  mcchain[:p] = Dict{Int, Matrix{Float64}}()
  for i in 1:nsites
    mcchain[:p][i] = Array(Float64, nmc, ncells[i])
  end

  counter::Int = 1
  tic()

  for k in 1:job[:gibbs].nsteps
    println("Running $k Gibbs iteration out of $(job[:gibbs].nsteps)...")
    tic()

    for i in 1:nsites
      wchain = Lora.run(job[:w][i](mcsample[:p][i], [mcsample[:w][i]]))
      mcsample[:w][i] = wchain.samples[end, :][1]
      @printf("  RAM acceptance rate for site %02d: %.2f%%\n", i, round(acceptance(wchain), 2))

      for j in 1:ncells[i]
        mcsample[:p][i][j] = run(job[:p][i][j](mcsample[:w][i]))
      end
    end

    if in(k, job[:gibbs].r)
      mcchain[:w][counter, :] = mcsample[:w]

      for i in 1:nsites
        for j in 1:ncells[i]
          mcchain[:p][i][counter, j] = mcsample[:p][i][j]
        end
      end

      counter += 1
    end

    println("  Gibbs step completed in $(toq()) seconds")
  end

  println("Total ellpased time: $(toq()) seconds")

  return mcchain
end
