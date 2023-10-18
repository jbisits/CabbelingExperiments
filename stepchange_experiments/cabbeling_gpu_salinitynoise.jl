using TwoLayerDirectNumericalShenanigans

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))
reduced_domain = (Lx = 0.1, Ly = 0.1, Lz = 0.8)
interface_location = 0.3
resolution = (Nx = 124, Ny = 124, Nz = 1100)

## Setup the dns_model
@info "Model setup"
dns_model = DNS(architecture, reduced_domain, resolution, diffusivities;
                reference_density = REFERENCE_DENSITY, zgrid_stretching = false)

## set initial conditions
@info "Setting initial conditions"
T₀ᵘ = -1.5
S₀ᵘ = 34.58
cabbeling = CabbelingUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
depth = find_depth(dns_model, interface_location)
profile_function = StepChange(depth)

## Salinity noise
depths = find_depth(dns_model, [interface_location+ 0.02, interface_location - 0.02])
scales = similar(depths)
fill!(scales, 2e-4)
initial_noise = SalinityNoise(depths, scales)
@info "Building DNS"
dns = TwoLayerDNS(dns_model, profile_function, initial_conditions; initial_noise)

@info "Setting two layer initial conditions"
set_two_layer_initial_conditions!(dns)

## build the simulation
Δt = 1e-4
stop_time = 3 * 60 * 60 # seconds
save_schedule = 60  # seconds
simulation = DNS_simulation_setup(dns, Δt, stop_time, save_schedule)

## Run the simulation
run!(simulation)
