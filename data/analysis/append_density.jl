using DirectNumericalCabbelingShenanigans.OutputUtilities

saved_simulations = readdir(SIMULATION_PATH, join = true)
for simulation ∈ saved_simulations
    open_sim = jldopen(simulation)
    if "σ₀" ∉ keys(open_sim["timeseries"])
        close(open_sim)
        compute_density!(simulation, density_string = "σ₀")
    else
        @info "A timeseries with that symbol already exists in $simulation."
        close(open_sim)
    end
end
