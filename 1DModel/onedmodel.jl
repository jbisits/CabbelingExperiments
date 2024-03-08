"""
    module OneDModel
One dimensional model of water column at lower resolution to compare DNS with. This model
has a convective adjustment scheme so we can compare the onset of instability with the DNS.
"""
module OneDModel

export run_OneDModel, hovmoller_plot, TShovmoller_plot

using Oceananigans, SeawaterPolynomials.TEOS10, GibbsSeaWater, CairoMakie
using Oceananigans.Units: seconds, minutes, hours, days
using Oceananigans: ∂z_b
using Oceananigans: Models.seawater_density

"""
    function run_OneDModel(salinity_initial_condition::Symbol; params)
Run a 1D model with initial conditions that are set to cabbeling instability at interface
between upper (quasi mixed, ``z ∈ (-100, 0)``) and lower region
(increasing salinity gradient, ``z ∈ (-500, -100)``).
There is a very slight temperature gradient in upper region to avoid spurios mixing.
By default, salinity and temperature have equal diffusivity but this can be changed by
passing a keyword argument for `convective_κz` and `background_κz` that is a `NamedTuple`
for the `T` and `S` diffusivities.

Arguments:
- `salinity_initial_condition`, `Symbol` either `:stable`, `:cabbeling`, `:unstable` for the
salinity in the upper layer.
"""
function run_OneDModel(salinity_initial_condition::Symbol;
                    Tᵤ = -1.5,
                    Tₗ = 0.5,
                    Sₗ = 34.7,
                    Sₘ = 34.7,   # maximum salinity
                    Nz = 1400,     # number of levels
                    Lz = -1000,    # overall depth
     reference_density = gsw_rho(Sₗ, Tₗ, 0),
         convective_κz = 10.0,
         background_κz = 1e-7,
                     ν = 1e-6,
   reference_gp_height = 0,
                    Δt = 1,       # minutes
              savepath = "OneDModelOutput",
            sim_length = 11 / 24,       # in days
             save_freq = 1        # in minutes
        )

    # Grid
    grid = RectilinearGrid(size = Nz, z = (Lz, 0), topology=(Flat, Flat, Bounded))

    # Buoyancy, using TEOS10
    EOS = TEOS10EquationOfState(; reference_density)
    buoyancy = SeawaterBuoyancy(equation_of_state = EOS)

    ## Turbulence closure
    closure = ConvectiveAdjustmentVerticalDiffusivity(background_νz = ν, convective_νz = ν;
                                                       convective_κz, background_κz)

    # Set temperature initial condition
    T₀ = Array{Float64}(undef, size(grid))
    # This adds a temperature gradient to avoid spurios convective mixing in the mixed layer
    Tₗ_array = fill(Tₗ, Int(Nz / 2))
    Tᵤ_array = reverse(range(Tᵤ, Tᵤ, length = Int(Nz / 2)))
    T₀[:, :, :] = vcat(Tₗ_array, Tᵤ_array)

    # Set the salinity initial condition
    S₀ = Array{Float64}(undef, size(grid))
    Sᵤ = getfield(salinity_initial_conditions, salinity_initial_condition)
    Sᵤ_array = fill(Sᵤ, Int(Nz / 2))
    Sₗ_array = range(Sₘ, Sₗ, length = Int(Nz / 2))
    S₀[:, :, :] = vcat(Sₗ_array, Sᵤ_array)

    savefile = savepath*"_"*string(salinity_initial_condition)*".jld2"

    @info "Setting up model and building simulation."
    model = NonhydrostaticModel(grid = grid,
                                tracers = (:T, :S),
                                buoyancy = buoyancy,
                                closure = closure)

    σ = seawater_density(model, geopotential_height = reference_gp_height)
    set!(model, T = T₀, S = S₀)

    simulation = Simulation(model, Δt = Δt * minutes, stop_time = sim_length * days)

    outputs = (T = model.tracers.T, S = model.tracers.S,
               κ = save_diffusivity, σ = σ)

    simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                            filename = savefile,
                                            schedule = TimeInterval(save_freq * minutes))

    run!(simulation)

    return nothing

end
"""
    const salinity_initial_conditions
Salinity initial conditions for upper layer temperature of `Tᵤ = -1.5°C`.
"""
const salinity_initial_conditions = (stable = 34.551, cabbeling = 34.58, unstable = 34.59,
                                     isohaline = 34.7, isothermal = 34.69431424)
"""
    save_∂z_b(model)
Save the buoyancy gradient that is calculated during the simulations.
The `ConvectiveAdjustmentVerticalDiffusivity` closure computes the buoyancy gradient using
`∂z_b` and applies the convective diffusivity if `∂z_b < 0` and the background diffusivity
if `∂z_b ≥ 0`.
"""
function save_∂z_b(model)

    buoyancy_gradient = Array{Float64}(undef, model.grid.Nz)

    for i ∈ 1:model.grid.Nz
        buoyancy_gradient[i] = ∂z_b(1, 1, i, model.grid, model.buoyancy, model.tracers)
    end

    return buoyancy_gradient
end
"""
    save_diffusivity(model)
Save the diffusivity that is applied to a region using the method in the
`ConvectiveAdjustmentVerticalDiffusivity` closure. The code obtained to do this is from
`Oceananigans.jl` though I do not use the `@kernel` macro as I am not on a `GPU`.
"""
function save_diffusivity(model)

    is_stableᶜᶜᶠ(i, j, k, grid, tracers, buoyancy)=∂z_b(i, j, k, grid, buoyancy, tracers)>=0
    diffusivities = Array{Float64}(undef, model.grid.Nz)

    for i ∈ 1:model.grid.Nz
        stable_cell = is_stableᶜᶜᶠ(1, 1, i, model.grid, model.tracers, model.buoyancy)

        diffusivities[i] = ifelse(stable_cell,
                                  model.closure.background_κz,
                                  model.closure.convective_κz)
    end

    return diffusivities
end
"""
    function hovmoller_plot(field::FieldTimeSeries, fieldname::AbstractString;
                            zrange = nothing,
                            colormap = :thermal)
Hovmoller plot of a `field`. Can zoom in region by passing `zrange` argument,
which is index of where you want to plot.
"""
function hovmoller_plot(field::FieldTimeSeries, fieldname::AbstractString;
                        zrange = nothing,
                        colormap = :thermal)

    z = isnothing(zrange) ? znodes(field[1]) : znodes(field[1])[zrange]
    t = field.times / (60 * 60) # hours
    plot_field = isnothing(zrange) ? interior(field, 1, 1, :, :)' :
                                     interior(field, 1, 1, zrange, :)'

    fig = Figure(size = (1000, 600))
    ax = Axis(fig[1, 1],
            xlabel = "t (hours)",
            xaxisposition = :top,
            ylabel = "z (m)")
    hm_s = heatmap!(ax, t, z, plot_field; colormap)
    Colorbar(fig[2, 1], hm_s, label = fieldname, vertical = false, flipaxis = false)

    return fig

end
"""
    function TShovmoller_plot(S_ts::FieldTimeSeries, T_ts::FieldTimeSeries; zrange = nothing)
Hovmoller plot of temperature and salinity from saved ouput of the 1D model. Can zoom in
region by passing `zrange` argument, which is index of where you want to plot.
"""
function TShovmoller_plot(S_ts::FieldTimeSeries, T_ts::FieldTimeSeries; zrange = nothing)

    z = isnothing(zrange) ? znodes(S_ts[1]) : znodes(S_ts[1])[zrange]
    t = S_ts.times / (60 * 60) # hours
    S = isnothing(zrange) ? interior(S_ts, 1, 1, :, :)' : interior(S_ts, 1, 1, zrange, :)'
    T = isnothing(zrange) ? interior(T_ts, 1, 1, :, :)' : interior(T_ts, 1, 1, zrange, :)'

    fig = Figure(size = (1000, 600))
    ax = [Axis(fig[1, i],
            xlabel = "t (hours)",
            xaxisposition = :top,
            ylabel = "z (m)") for i ∈ 1:2]
    linkyaxes!(ax[1], ax[2])
    hm_s = heatmap!(ax[1], t, z, S; colormap = :haline)
    Colorbar(fig[2, 1], hm_s, label = "Salinity (g/kg)", vertical = false, flipaxis = false)
    hm_t = heatmap!(ax[2], t, z, T; colormap = :thermal)
    Colorbar(fig[2, 2], hm_t, label = "Temperature (°C)", vertical = false, flipaxis = false)

    return fig
end

end # module
