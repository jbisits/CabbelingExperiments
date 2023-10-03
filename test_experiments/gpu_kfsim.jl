## A test model on the GPU
architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))
DNS_resolution = (Nx = 100, Ny = 100, Nz = 1000)

## Setup the model
@info "Model setup"
model = DNS(architecture, DOMAIN_EXTENT, DNS_resolution, diffusivities;
            reference_density = REFERENCE_DENSITY, zgrid_stretching = false)

## set initial conditions
@info "Setting initial conditions"
T₀ᵘ = -1.5
S₀ᵘ = 34.58
cabbeling = CabbelingUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
depth = find_depth(model, INTERFACE_LOCATION)
profile_function = StepChange(depth)

## Salinity noise
depths = find_depth(model, [INTERFACE_LOCATION + 0.02, INTERFACE_LOCATION - 0.02])
scales = similar(depths)
fill!(scales, 2e-4)
initial_noise = SalinityNoise(depths, scales)
@info "Building DNS"
dns = TwoLayerDNS(model, profile_function, initial_conditions; initial_noise)

@info "Setting two layer initial conditions"
set_two_layer_initial_conditions!(dns)

## build the simulation
Δt = 1e-4
stop_time = 3 * 60
save_schedule = 1 # seconds
simulation = DNS_simulation_setup_test(dns, Δt, stop_time, save_schedule)

## Run the simulation
run!(simulation)
