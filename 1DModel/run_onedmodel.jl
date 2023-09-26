include("onedmodel.jl")
using .OneDModel

run_OneDModel(:cabbeling)

##
using TwoLayerDirectNumericalShenanigans, CairoMakie

## Output
output = joinpath(@__DIR__, "OneDModelOutput_cabbeling.jld2")
# compute_density!(output)
S_ts, T_ts = FieldTimeSeries(output, "S"), FieldTimeSeries(output, "T")
σ₀_ts =  FieldTimeSeries(output, "σ")
w = FieldTimeSeries(output, "w")
κ_interface = [jldopen(output)["timeseries"]["κ"]["$t"][400] for t ∈ string.(0:1440)]

## Initial snapshots
visualise_snapshot(T_ts, "Θ (°C)", 1, 1, 1)
visualise_snapshot(S_ts, "S (gkg⁻ꜝ)", 1, 1, 1; colormap = :haline)

## Hovmoller plots
hplot = TShovmoller_plot(S_ts, T_ts; zrange = 300:500)
hplot_density = hovmoller_plot(σ₀_ts, "σ₀ (kgm⁻³)"; colormap = :dense, zrange = 400:500)

## Diffusivity
lines(t, κ_interface)
findfirst(κ_interface .== 1)
t[10] # so convective adjustment scheme triggered after 0.15 hours ⟹ after 9 mins, 540s
