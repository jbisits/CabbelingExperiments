# High resolution two layer simulation
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

architecture = GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-6, T = 1e-5))
resolution = (Nx = 20, Ny = 20, Nz = 4000)

## Setup the model
model = DNS(architecture, domain_extent, resolution, diffusivities; reference_density)

## set initial conditions
T₀ᵘ = -1.5
S₀ᵘ = 34.551
stable = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(stable)
set_two_layer_initial_conditions!(model, initial_conditions, interface_location,
                                  salinity_perturbation = true, t = 0.0001)

## build the simulation
Δt = 1e-5
stop_time = 2
save_schedule = 0.2 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
