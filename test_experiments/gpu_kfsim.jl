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
    parameters = (pᵣ = density_reference_pressure,)
    σ = DensityField(model, parameters)

    # Inferred vertical diffusivity
    b_anomaly = BuoyancyField(model)
    b_zgrad = Field(∂b∂z(model))

    wᶜᶜᶜ_anomaly = Field(wᶜᶜᶜ(dns.model))
    κᵥ = Integral((-wᶜᶜᶜ_anomaly * b_anomaly) / b_zgrad)

    # Minimum in space Kolmogorov length scale
    ϵ = KineticEnergyDissipationRate(model)
    ∫ϵ = Integral(ϵ)
    η_space(model) = (model.closure.ν^3 / maximum(ϵ))^(1/4)

    # Dimensions and attributes for custom saved output
    dims = Dict("η_space" => (), "σ" => ("xC", "xC", "zC"), "κᵥ" => (), "∫ϵ" => ())
    oa = Dict(
        "σ" => Dict("longname" => "Seawater potential density calculated using TEOS-10 at $(density_reference_pressure)dbar",
                    "units" => "kgm⁻³"),
        "η_space" => Dict("longname" => "Minimum (in space) Kolmogorov length"),
        "κᵥ" => Dict("longname" => "Inferred vertical diffusivity",
                     "units" => "m²s⁻¹"),
        "∫ϵ" => Dict("longname" => "Volume integrated turbulent kintetic energy dissipation")
        )

    # outputs to be saved during the simulation
    outputs = Dict("S" => S, "T" => T, "η_space" => η_space, "σ" => σ, "κᵥ" => κᵥ, "∫ϵ" => ∫ϵ)
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
simulation_progress(sim) = @printf("i: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
            iteration(sim), time(sim), prettytime(sim.run_wall_time),
            sim.Δt, AdvectiveCFL(sim.Δt)(sim.model),
            DiffusiveCFL(sim.Δt)(sim.model))

## build the simulation
Δt = 1e-4
stop_time = 3 * 60
save_schedule = 1 # seconds
simulation = DNS_simulation_setup_test(dns, Δt, stop_time, save_schedule)

## Run the simulation
run!(simulation)
