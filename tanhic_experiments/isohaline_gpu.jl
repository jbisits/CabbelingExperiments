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
isohaline = IsohalineUpperLayerInitialConditions(T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(isohaline)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 100.0)
set_two_layer_initial_conditions!(model, initial_conditions, profile_function)

## build the simulation
Δt = 1e-5
stop_time = 5
save_schedule = 1 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
