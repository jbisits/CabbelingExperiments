include("onedmodel.jl")
using .OneDModel

run_OneDModel(:isothermal, background_κz = 1e-2, Tᵤ = 0.5, salinity_noise = true, Sgrad = 2e-5)
