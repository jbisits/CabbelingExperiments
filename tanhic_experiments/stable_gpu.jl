# High resolution two layer simulation
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

architecture = GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions
T₀ᵘ = -1.5
S₀ᵘ = 34.551
stable = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(stable)
interface_width = 75
set_two_layer_initial_conditions!(model, initial_conditions, INTERFACE_LOCATION, :tanh,
                                  interface_width;
                                  salinity_perturbation = true,
                                  salinity_perturbation_width = 30)
add_velocity_random_noise!(model, 1e-2, INTERFACE_LOCATION / 2)

## build the simulation
Δt = 1e-5
stop_time = 10 # seconds (in simulation time)
save_schedule = 1 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
