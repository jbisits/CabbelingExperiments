using DirectNumericalCabbelingShenanigans
using Oceananigans.Fields
using DirectNumericalCabbelingShenanigans.OutputUtilities

## Load in saved output
# sim_path = joinpath(SIMULATION_PATH,
#                     "tanhic/equal_diffusivity/perturb_salinity/cabbeling.jld2") # local
sim_path = joinpath(SIMULATION_PATH, "cabbeling.jld2") # gadi
T_ts = FieldTimeSeries(sim_path, "T")
S_ts = FieldTimeSeries(sim_path, "S")

## Initial snapshots
visualise_snapshot(T_ts, "Θ (°C)", 10, 10, 11)
visualise_snapshot(S_ts, "S (gkg⁻ꜝ)", 10, 10, 11; colormap = :haline)

## Animations (x-z)
animate_2D_field(T_ts, "Θ (°C)", 1, 1)

## Compute a density `FieldTimeSeries`
σ₀_ts = compute_density(S_ts, T_ts)
visualise_snapshot(σ₀_ts, "σ₀ (kgm⁻³)", 10, 10, 10; colormap = :dense)
animate_2D_field(σ₀_ts, "σ₀ (kgm⁻³)", 10, 10; colormap = :dense)
