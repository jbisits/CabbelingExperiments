using TwoLayerDirectNumericalShenanigans

restart = true

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))
domain_extent = (Lx = 0.07, Ly = 0.07, Lz = -1.0)
resolution = (Nx = 35, Ny = 35, Nz = 500)
eos = TEOS10EquationOfState(reference_density = REFERENCE_DENSITY)

## Setup the dns_model
@info "Model setup"
dns_model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

## set initial conditions
@info "Setting initial conditions"
T₀ᵘ = 0.5
S₀ᵘ = 34.69431424
const interface_location = -0.5
isothermal = IsothermalUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(isothermal)
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
Δt = 1e-2
max_Δt = 7e-2
stop_time = 13 * 60 * 60 # seconds
save_schedule = 60  # seconds
checkpointer_time_interval = 60 * 60 # seconds
output_path = @__DIR__
@info "Setting up simulation"
simulation = TLDNS_simulation_setup(tldns, Δt, stop_time, save_schedule, TLDNS.save_computed_output!,
                                    TLDNS.save_vertical_velocities!;
                                    checkpointer_time_interval, output_path, max_Δt,
                                    overwrite_saved_output = restart)

pickup = restart ? false : true
## Run the simulation
run!(simulation; pickup)
