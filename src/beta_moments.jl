beta_a_from_mv(m::Float64, v::Float64) = (m*(1-m)/v-1)*m

beta_b_from_mv(m::Float64, v::Float64) = (m*(1-m)/v-1)*(1-m)

function beta_pars_from_mv(m::Float64, v::Float64)
  c = m*(1-m)/v-1
  return c*m, c*(1-m)
end
