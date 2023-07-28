# High resolution two layer simulation
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

architecture = GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-6, T = 1e-5))

## Setup the model
model = DNS(architecture, domain_extent, high_resolution, diffusivities; reference_density)

## set initial conditions
T₀ᵘ = -1.5
isohaline = IsohalineUpperLayerInitialConditions(T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(isohaline)
start_time = 0.1
set_two_layer_initial_conditions!(model, initial_conditions, interface_location, start_time)

## build the simulation
Δt = 1e-5
stop_time = 2
save_schedule = 1 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
