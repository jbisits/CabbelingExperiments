using Oceananigans.Fields

## Load (lazily) in saved output
sim_path = joinpath(SIMULATION_PATH, "verystable_tanh_salinitygaussianprofile_mp_10min.jld2")
T_ts = FieldTimeSeries(sim_path, "T", backend = OnDisk())
S_ts = FieldTimeSeries(sim_path, "S", backend = OnDisk())
# σ₀_ts = FieldTimeSeries(sim_path, "σ₀", backend = OnDisk())
