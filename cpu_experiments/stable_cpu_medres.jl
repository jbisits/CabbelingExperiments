# High resolution two layer simulation
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

architecture = CPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))
resolution = (Nx = 50, Ny = 50, Nz = 1400)

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, resolution, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions
T₀ᵘ = -1.5
S₀ᵘ = 34.551
stable = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(stable)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 1000.0)
z = znodes(model.grid, Center(), Center(), Center())
depth_idx = findfirst(z .> 0.98 * INTERFACE_LOCATION)
salinity_perturbation = GaussianBlob(z[depth_idx], [0.0, 0.0], 0.85)
set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
                                  salinity_perturbation)
DNCS.OutputUtilities.visualise_initial_conditions(model, 1, 1)
DNCS.OutputUtilities.visualise_initial_density(model, 1, 1, 0)
## build the simulation
Δt = 1e-5
stop_time = 5 * 60 # seconds (in simulation time)
save_schedule = 5 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions, max_Δt = 0.2)

## Run the simulation
run!(simulation)
