using DirectNumericalCabbelingShenanigans
using Oceananigans.Fields
using DirectNumericalCabbelingShenanigans.OutputUtilities

## Load in saved output
sim_path = joinpath(SIMULATION_PATH, "tanhic/equal_diffusivity/stable_medres_1_5min.jld2")
T_ts = FieldTimeSeries(sim_path, "T")
S_ts = FieldTimeSeries(sim_path, "S")

## Initial snapshots
visualise_snapshot(T_ts, "Θ (°C)", 1)
visualise_snapshot(S_ts, "S (gkg⁻ꜝ)", 1; colormap = :haline)

## Animations (x-z)
animate_2D_field(T_ts, "Θ (°C)", (:x, :z))

## Compute a density `FieldTimeSeries`
σ₀_ts = compute_density(S_ts, T_ts)
visualise_snapshot(σ₀_ts, "σ₀ (kgm⁻³)", 1; colormap = :dense)
animate_2D_field(σ₀_ts, "σ₀ (kgm⁻³)", (:x, :z); colormap = :dense)

## temp till I fix `animate_2D_field`
x, y, z = nodes(σ₀_ts[1])
t = σ₀_ts.times
n = Observable(1)
field_tₙ = @lift interior(σ₀_ts[$n], :, 1, :)
c_limits = extrema(interior(σ₀_ts, :, :, :, 1)) .+ (-0.001, 0.001)
title = @lift @sprintf("t=%1.2f", t[$n])
fig, ax, hm = heatmap(x, z, field_tₙ;
                      colormap = :dense, colorrange = c_limits)
ax.aspect = 0.25
ax.xlabel = "x"
ax.ylabel = "z"
Colorbar(fig[1, 2], hm, label = "σ₀ (kgm⁻³)")

frames = eachindex(t)
filename =  "xy_σ₀ (kgm⁻³)"
record(fig, joinpath(pwd(), filename * ".mp4"),
      frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")
    n[] = i
end
