### This is not a generically coded Metropolis-within-Gibbs sampler, it is tailored to the problem at hand

function metropolis_within_gibbs(m::Vector{Float64}, n::Matrix{Float64}, x::Matrix{Float64})
  # Indices (i, j) refer to the i-th site (transcript) on the j-th cell
  # data_i = (mean_i, {n_ij: j}, {x_ij: j}), mean_i = x_i/n_i
  # Step 1: ∀ i draw a sample from v_i | ({p_ij: j}, data_i) using Metropolis sampling
  # Step 2: ∀ (i, j) draw a sample from p_ij | (v_i, data_i) which follows a Beta distribution

  nsites, ncells = size(n)
  @assert length(m) == nsites "length(m) must be equal to size(n, 1)"
  @assert size(x) == (nsites, ncells) "size(x) must be equal to size(n)"

  v = Array(Float64, nsites)
  p = Array(Float64, nsites, ncells)

  for i in 1:nsites
    v[i] = metropolis()

    map(sample_from_conjugate_beta(), )
  end
end
