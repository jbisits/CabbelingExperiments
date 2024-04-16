include("onedmodel.jl")
using .OneDModel

run_OneDModel(:isothermal; Tᵤ = 0.5, Nz = 100, background_κz = 1e-5)
