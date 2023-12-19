using TwoLayerDirectNumericalShenanigans

restart = true

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
eos = TEOS10EquationOfState(reference_density = REFERENCE_DENSITY)

## Setup the dns_model
@info "Model setup"
dns_model = DNSModel(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities, eos)

## set initial conditions
@info "Setting initial conditions"
T₀ᵘ = -1.5
S₀ᵘ = 34.551
stable = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(stable)
depth = find_depth(dns_model, INTERFACE_LOCATION)
profile_function = StepChange(depth)

## Salinity noise
depths = find_depth(dns_model, [INTERFACE_LOCATION + 0.02, INTERFACE_LOCATION - 0.02])
scales = similar(depths)
fill!(scales, 2e-4)
initial_noise = SalinityNoise(depths, scales)
@info "Building DNS"
tldns = TwoLayerDNS(dns_model, profile_function, initial_conditions; initial_noise)

@info "Setting two layer initial conditions"
set_two_layer_initial_conditions!(tldns)

## build the simulation
Δt = 1e-4
max_Δt = 0.1
stop_time = 60 * 60 # seconds
save_schedule = 60  # seconds
checkpointer_time_interval = 30 * 60 # seconds
output_path = joinpath(@__DIR__, "outputs_doublediffusion/")
@info "Setting up simulation"
simulation = TLDNS_simulation_setup(tldns, Δt, stop_time, save_schedule, TLDNS.save_computed_output!;
                                    checkpointer_time_interval, output_path, max_Δt,
                                    overwrite_saved_output = restart)

pickup = restart ? false : true
## Run the simulation
run!(simulation; pickup)
