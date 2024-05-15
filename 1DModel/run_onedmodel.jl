include("onedmodel.jl")
using .OneDModel

run_OneDModel(:isothermal, Tᵤ = 0.5, salinity_noise = true)
run_OneDModel(:cabbeling, Tᵤ = -1.5, convective_κz = 1e-4, salinity_noise = true)
