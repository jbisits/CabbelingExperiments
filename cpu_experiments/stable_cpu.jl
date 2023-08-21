# High resolution two layer simulation
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

architecture = CPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions
T₀ᵘ = -0.5
S₀ᵘ = 34.625
stable = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(stable)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 3500.0)
## `GaussianBlob`
z = znodes(model.grid, Center(), Center(), Center())
depth_idx = findfirst(z .> 9 * INTERFACE_LOCATION / 10)
salinity_perturbation = GaussianBlob(z[depth_idx], [0.0, 0.0], 10.0)
set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
                                  salinity_perturbation)
## `GaussianProfile` with `RandomPerturbations`
salinity_perturbation = GaussianProfile(INTERFACE_LOCATION, INTERFACE_LOCATION / 1.1,
                                        100.0, 2.0)
z = znodes(model.grid, Center(), Center(), Center())
depth_idx = findfirst(z .> INTERFACE_LOCATION / 1.1)
z[depth_idx]
salinity_noise = RandomPerturbations(-0.34077380952380953, 0.1)
set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
                                  salinity_perturbation, salinity_noise)
## Look at the output
DNCS.OutputUtilities.visualise_initial_conditions(model, 1, 1)
DNCS.OutputUtilities.visualise_initial_density(model, 1, 1, 0)

## build the simulation
Δt = 1e-5
stop_time = 60 # seconds (in simulation time)
save_schedule = 1 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
