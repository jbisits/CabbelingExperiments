using TwoLayerDirectNumericalShenanigans

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions
T₀ᵘ = -1.5
S₀ᵘ = 34.568
cabbeling = CabbelingUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
depth = find_depth(model, INTERFACE_LOCATION)
profile_function = StepChange(depth)

## Salinity noise
depths = find_depth(model, [INTERFACE_LOCATION + 0.02, INTERFACE_LOCATION - 0.02])
scales = similar(depths)
fill!(scales, 2e-4)
initial_noise = SalinityNoise(depths, scales)
dns = TwoLayerDNS(model, profile_function, initial_conditions; initial_noise)

set_two_layer_initial_conditions!(dns)

## build the simulation
Δt = 1e-4
stop_time = 10 * 60
save_schedule = 5 # seconds
simulation = DNS_simulation_setup(dns, Δt, stop_time, save_schedule)

## Run the simulation
run!(simulation)
