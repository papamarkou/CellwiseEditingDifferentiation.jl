function metropolis_within_gibbs(m::Float64, n::Vector{Float64}, x::Vector{Float64})
  # Indices (i, j) refer to the i-th site (transcript) on the j-th cell
  # Step 1: ∀ i draw a sample from v_i | {p_ij} using Metropolis sampling
  # Step 2: ∀ (i, j) draw a sample from p_ij | v_i, which follows a Beta distribution
end
