### This is not a generically coded Metropolis-within-Gibbs sampler, it is tailored to the problem at hand

function gibbs(data::Dict{Symbol, Any}, init::Dict{Symbol, Any}, job)
  nsites, ncells = size(data[:coverage])
  @assert length(data[:m]) == nsites "length(data[:m]) must be equal to size(data[:coverage], 1)"
  @assert size(data[:edited]) == (nsites, ncells) "size(data[:edited]) must be equal to size(data[:coverage])"

  nsamples = length(job[:gibbs].r)

  mcsample = Dict{Symbol, Any}(:w=>init[:w], :p=>init[:p])
  mcchain = Dict{Symbol, Any}(:w=>Array(Float64, nsamples, nsites), :p=>Array(Float64, nsamples, nsites, ncells))

  counter::Int = 1
  tic()

  for k in 1:job[:gibbs].nsteps
    print("Running $k Gibbs iteration out of $(job[:gibbs].nsteps)...")
    tic()

    for i in 1:nsites
      mcsample[:w][i] = Lora.run(job[:w][i](vec(mcsample[:p][i, :]), [mcsample[:w][i]])).samples[end, :][1]

      for j in 1:ncells
        mcsample[:p][i, j] = run(job[:p][i, j](mcsample[:w][i]))
      end
    end

    if in(k, job[:gibbs].r)
      mcchain[:w][counter, :] = mcsample[:w]
      mcchain[:p][counter, :, :] = mcsample[:p]
      counter += 1
    end

    println(" completed in $(toq()) seconds")
  end

  println("Total ellpased time: $(toq()) seconds")

  return mcchain
end
