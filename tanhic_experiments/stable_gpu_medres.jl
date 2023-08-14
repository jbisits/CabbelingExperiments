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
S₀ᵘ = 34.551
stable = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(stable)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 100.0)
z = znodes(model.grid, Center(), Center(), Center())
#depth_idx = findfirst(z .> 2 * INTERFACE_LOCATION / 3)
salinity_perturbation = GaussianBlob(-0.24958333333333338, [0.0, 0.0], 1.5)
set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
                                  salinity_perturbation)

## build the simulation
Δt = 1e-5
stop_time = 10 * 60 # seconds (in simulation time)
save_schedule = 1 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
