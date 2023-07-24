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

## Hovmoller plot
z = znodes(S_ts[1])
t = S_ts.times / (60 * 60) # hours

fig = Figure(size = (1000, 600))
ax = [Axis(fig[1, i],
           xlabel = "t (hours)",
           xaxisposition = :top,
           ylabel = "z (m)") for i ∈ 1:2]
linkyaxes!(ax[1], ax[2])
hm_s = heatmap!(ax[1], t, z, interior(S_ts, 1, 1, :, :)', colormap = :haline)
Colorbar(fig[2, 1], hm_s, label = "Salinity (g/kg)", vertical = false, flipaxis = false)
hm_t = heatmap!(ax[2], t, z, interior(T_ts, 1, 1, :, :)', colormap = :thermal)
Colorbar(fig[2, 2], hm_t, label = "Temperature (°C)", vertical = false, flipaxis = false)

fig

## Diffusivity
lines(t, κ_interface)
findfirst(κ_interface .== 1)
t[10] # so convective adjustment scheme triggered after 0.15 hours ⟹ after 9 mins, 540s
