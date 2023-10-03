function DNS_simulation_setup_test(dns::TwoLayerDNS, Δt::Number,
                            stop_time::Number, save_schedule::Number,
                            output_writer::Symbol=:netcdf;
                            cfl = 0.75,
                            diffusive_cfl = 0.75,
                            max_change = 1.2,
                            max_Δt = 1e-1,
                            density_reference_pressure = 0,
                            save_velocities = false)

    model = dns.model
    simulation = Simulation(model; Δt, stop_time)

    # time step adjustments
    wizard = TimeStepWizard(; cfl, diffusive_cfl, max_change, max_Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    # model tracers
    S, T = model.tracers.S, model.tracers.T
    # custom saved output

    # Density
    σ = TLDNS.DensityField(model, density_reference_pressure)

    # Inferred vertical diffusivity
    σ_anomaly_interpolated = TLDNS.InterpolatedDensityAnomaly(model, density_reference_pressure)
    w = model.velocities.w
    κᵥ = Integral((-w * σ_anomaly_interpolated) / σ)

    # Minimum in space Kolmogorov length scale
    #ϵ = KineticEnergyDissipationRate(model)
    #η_space(model) = minimum(model.closure.ν ./ ϵ)

    # Dimensions and attributes for custom saved output
    dims = Dict("σ" => ("xC", "xC", "zC"), "κᵥ" => ())
    oa = Dict(
                "σ" => Dict("longname" => "Seawater potential density calculated using TEOS-10 at $(density_reference_pressure)dbar",
                "units" => "kgm⁻³"),
                #"η_space" => Dict("longname" => "Minimum (in space) Kolmogorov length"),
                "κᵥ" => Dict("longname" => "Inferred vertical diffusivity",
                "units" => "m²s⁻¹"))

    # outputs to be saved during the simulation
    outputs = Dict("S" => S, "T" => T, "σ" => σ, "κᵥ" => κᵥ)
    if save_velocities
        u, v = model.velocities.u, model.velocities.v
        velocities = Dict("u" => u, "v" => v, "w" => w)
        merge!(outputs, velocities)
    end

    filename = TLDNS.form_filename(dns, stop_time, output_writer)
    simulation.output_writers[:outputs] = output_writer == :netcdf ?
                    NetCDFOutputWriter(model, outputs,
                                    filename = filename,
                                    schedule = TimeInterval(save_schedule),
                                    overwrite_existing = true,
                                    dimensions = dims,
                                    output_attributes = oa
                                    ) :
                    JLD2OutputWriter(model, outputs,
                                    filename = filename,
                                    schedule = TimeInterval(save_schedule),
                                    overwrite_existing = true)

    TLDNS.non_dimensional_numbers!(simulation, dns)
    TLDNS.predicted_maximum_density!(simulation, dns)

    # progress reporting
    simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(100))

    return simulation

end
