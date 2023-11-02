using TwoLayerDirectNumericalShenanigans

restart = true

architecture = GPU()
diffusivities = (ν = 5e-6, κ = (S = 5e-9, T = 5e-7))
eos = TEOS10EquationOfState(reference_density = REFERENCE_DENSITY)

## Setup the dns_model
@info "Model setup"
dns_model = DNSModel(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities, eos)

## set initial conditions
@info "Setting initial conditions"
T₀ᵘ = 0.5
S₀ᵘ = 34.69431424
isothermal = IsothermalUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(isothermal)
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
max_Δt = 0.075
stop_time = 20 * 60 # seconds
save_schedule = 30  # seconds
checkpointer_time_interval = 5 * 60 # seconds
output_path = joinpath(@__DIR__, "outputs_doublediffusion/")
@info "Setting up simulation"
simulation = TLDNS_simulation_setup(tldns, Δt, stop_time, save_schedule, TLDNS.save_computed_output!;
                                    checkpointer_time_interval, output_path, max_Δt,
                                    overwrite_saved_output = restart)

pickup = restart ? false : true
## Run the simulation
run!(simulation; pickup)
