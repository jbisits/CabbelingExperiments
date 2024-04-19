include("onedmodel.jl")
using .OneDModel

run_OneDModel(:isothermal; Tᵤ = 0.5, background_κz = 1e-2, ν = 0.0)
