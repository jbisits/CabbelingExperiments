include("onedmodel.jl")
using .OneDModel

run_OneDModel(:cabbeling, background_κz = 1e-2, ν = 0.0, salinity_noise = true)
