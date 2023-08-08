# High resolution two layer simulation
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))
resolution = (Nx = 20, Ny = 20, Nz = 1000)

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, resolution, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions
T₀ᵘ = -1.5
S₀ᵘ = 34.59
unstable = UnstableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(unstable)
interface_width = 100
set_two_layer_initial_conditions!(model, initial_conditions, INTERFACE_LOCATION, :tanh,
                                  interface_width)
add_velocity_random_noise!(model, 1e-2)

## build the simulation
Δt = 1e-5
stop_time = 60
save_schedule = 1 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
