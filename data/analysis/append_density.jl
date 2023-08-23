using DirectNumericalCabbelingShenanigans.OutputUtilities: compute_density!

saved_simulations = readdir(SIMULATION_PATH)
for simulation ∈ saved_simulations
    if "σ₀" ∉ keys(jldopen(simulation)["timeseries"])
        compute_density!(simulation, density_string = "σ₀")
    else
        @info "A timeseries with that symbol already exists in $simulation."
    end
end
