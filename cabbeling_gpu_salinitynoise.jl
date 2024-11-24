using TwoLayerDirectNumericalShenanigans
using Oceanostics: KineticEnergyDissipationRate, KineticEnergy
using Oceananigans: seawater_density

"Temporary to save the full epsilon field"
function save_computed_output_with_epsilon!(simulation, tldns, save_schedule, save_file, output_dir,
                                            overwrite_saved_output, reference_gp_height)

    model = tldns.model
    σ = seawater_density(model, geopotential_height = reference_gp_height)
    ϵ = KineticEnergyDissipationRate(model)
    ∫ϵ = Integral(ϵ)
    ϵ_maximum = Reduction(maximum!, ϵ, dims = (1, 2, 3))
    Eₖ = KineticEnergy(model)
    ∫Eₖ = Integral(Eₖ)

    computed_outputs = Dict("σ" => σ, "∫ϵ" => ∫ϵ, "ϵ_maximum" => ϵ_maximum, "∫Eₖ" => ∫Eₖ,
                            "∫Eₚ" => ∫Eₚ, "ϵ" => ϵ)

    oa = Dict(
        "σ" => Dict("longname" => "Seawater potential density calculated using TEOS-10 at $(reference_gp_height)dbar",
        "units" => "kgm⁻³"),
        "ϵ_maximum" => Dict("longname" => "Maximum (in space) TKE dissipation"),
        "∫ϵ" => Dict("longname" => "Volume integrated turbulent kintetic energy dissipation"),
        "∫Eₖ" => Dict("longname" => "Volume integrated turbulent kinetic energy"),
        "ϵ" => Dict("longname" => "TKE dissipation field")
    )
    simulation.output_writers[:computed_output] =
    save_file == :netcdf ? NetCDFOutputWriter(model, computed_outputs;
                        filename = "computed_output",
                        dir = output_dir,
                        overwrite_existing = overwrite_saved_output,
                        schedule = TimeInterval(save_schedule),
                        output_attributes = oa
                        ) :
        JLD2OutputWriter(model, computed_outputs;
                        filename = "computed_output",
                        dir = output_dir,
                        schedule = TimeInterval(save_schedule),
                        overwrite_existing = overwrite_saved_output)

    return nothing

end
save_computed_output_with_epsilon!(simulation, tldns, save_info::Tuple; reference_gp_height = 0) =
    save_computed_output_with_epsilon!(simulation, tldns, save_info..., reference_gp_height)


restart = true

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))
domain_extent = (Lx = 0.07, Ly = 0.07, Lz = -1.0)
resolution = (Nx = 115, Ny = 115, Nz = 1650)
eos = TEOS10EquationOfState(reference_density = REFERENCE_DENSITY)

## Setup the dns_model
@info "Model setup"
dns_model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

## set initial conditions
@info "Setting initial conditions in upper layer"
T₀ᵘ = -1.5
S₀ᵘ = 34.58
const interface_location = -0.5
cabbeling = CabbelingUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
depth = find_depth(dns_model, interface_location)
profile_function = StepChange(depth)

## Salinity noise
depths = find_depth(dns_model, [interface_location + 0.02, interface_location - 0.02])
scales = similar(depths)
fill!(scales, 2e-4)
initial_noise = SalinityNoise(depths, scales)
@info "Building DNS"
tldns = TwoLayerDNS(dns_model, profile_function, initial_conditions; initial_noise)

@info "Setting two layer initial conditions"
set_two_layer_initial_conditions!(tldns)

## build the simulation
Δt = 1e-3
max_Δt = 8e-2
stop_time = 3 * 60 * 60 # seconds
save_schedule = 60  # seconds
checkpointer_time_interval = 60 * 60 # seconds
output_path = joinpath(@__DIR__, "outputs_equaldiffusion/")
@info "Setting up simulation"

simulation = TLDNS_simulation_setup(tldns, Δt, stop_time, save_schedule, save_computed_output_with_epsilon!,
                                    TLDNS.save_vertical_velocities!;
                                    checkpointer_time_interval, output_path, max_Δt,
                                    overwrite_saved_output = restart)
pickup = restart ? false : true
## Run the simulation
run!(simulation; pickup)
