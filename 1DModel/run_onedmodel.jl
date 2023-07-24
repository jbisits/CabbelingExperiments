include("onedmodel.jl")
using .OneDModel

run_OneDModel(:cabbeling)
run_OneDModel(:cabbeling;
              ν = (S = 1e-4, T = 1e-4),
              background_κz = (S = 1e-6, T = 1e-5),
              convective_κz = (S = 1, T = 1))

## Output
output = joinpath(@__DIR__, "OneDModelOutput_cabbeling.jld2")
S_ts, T_ts = FieldTimeSeries(output, "S"), FieldTimeSeries(output, "T")
κ_interface = [jldopen(output)["timeseries"]["κ"]["$t"][400] for t ∈ string.(0:1440)]

using DirectNumericalCabbelingShenanigans.OutputUtilities

## Initial snapshots
visualise_snapshot(T_ts, "Θ (°C)", 1)
visualise_snapshot(S_ts, "S (gkg⁻ꜝ)", 1; colormap = :haline)

## Hovmoller plots
hplot = TShovmoller_plot(S_ts, T_ts; zrange = 300:500)
σ₀_ts = compute_density(S_ts, T_ts)
hplot_density = hovmoller_plot(σ₀_ts, "σ₀ (kgm⁻³)"; colormap = :dense, zrange = 300:500)

## Diffusivity
lines(t, κ_interface)
findfirst(κ_interface .== 1)
t[10] # so convective adjustment scheme triggered after 0.15 hours ⟹ after 9 mins, 540s
