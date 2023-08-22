# High resolution two layer simulation
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

architecture = GPU()
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
# depth_idx = findfirst(z .> INTERFACE_LOCATION / 2)
salinity_perturbation = GaussianBlob(-0.3560416666666667, [0.0, 0.0], 10.0)
set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
                                  salinity_perturbation)

## build the simulation
Δt = 1e-5
stop_time = 8 * 60 # seconds (in simulation time)
save_schedule = 5 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
