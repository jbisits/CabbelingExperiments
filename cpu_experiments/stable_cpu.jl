using TwoLayerDirectNumericalShenanigans

architecture = CPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions
T₀ᵘ = -1.5
S₀ᵘ = 34.551
stable = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(stable)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 3500.0)
## `GaussianBlob`
# z = znodes(model.grid, Center(), Center(), Center())
# depth_idx = findfirst(z .> INTERFACE_LOCATION / 1.1)
# depth = find_depth(model, INTERFACE_LOCATION / 1.1)
# tracer_perturbation = SalinityGaussianBlob(depth, [0.0, 0.0], 10.0)
# dns = TwoLayerDNS(model, profile_function, initial_conditions; tracer_perturbation)
# set_two_layer_initial_conditions!(dns)
## `GaussianProfile`
tracer_perturbation = SalinityGaussianProfile(INTERFACE_LOCATION, INTERFACE_LOCATION / 1.1,
                                                100.0, 8.0)

# set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
#                                   tracer_perturbation)

## With `RandomPerturbations`
# z = znodes(model.grid, Center(), Center(), Center())
# depth_idx = findfirst(z .> INTERFACE_LOCATION / 1.1)
depth = find_depth(model, INTERFACE_LOCATION / 1.1)
initial_noise = SalinityNoise(depth, 0.001)
dns = TwoLayerDNS(model, profile_function, initial_conditions;
                 tracer_perturbation, initial_noise)
set_two_layer_initial_conditions!(dns)

## Look at the output
using CairoMakie
visualise_initial_conditions(dns, 1, 1)
visualise_initial_density(dns, 1, 1, 0)

## build the simulation
Δt = 1e-5
stop_time = 60 # seconds (in simulation time)
save_schedule = 1 # seconds
simulation = DNS_simulation_setup(dns, Δt, stop_time, save_schedule)

## Run the simulation
run!(simulation)
