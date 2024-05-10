include("onedmodel.jl")
using .OneDModel

run_OneDModel(:isothermal, Tᵤ = 0.5, salinity_noise = true, Sgrad = 2e-5)
