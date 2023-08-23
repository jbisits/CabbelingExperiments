using Oceananigans.Fields
using DirectNumericalCabbelingShenanigans.OutputUtilities: compute_density

## Load in saved output
sim_path = joinpath(SIMULATION_PATH, "stable.jld2")
T_ts = FieldTimeSeries(sim_path, "T", backend = OnDisk())
S_ts = FieldTimeSeries(sim_path, "S", backend = OnDisk())

## Compute a density `FieldTimeSeries`
σ₀_ts = compute_density(S_ts, T_ts)
