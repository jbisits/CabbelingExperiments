# High resolution two layer simulation
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))
resolution = (Nx = 50, Ny = 50, Nz = 1400)

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, resolution, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions
T₀ᵘ = -1.5
S₀ᵘ = 34.568
cabbeling = CabbelingUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 1000.0)
# Need to find this depth manually because `GPU()` does not like scalar indexing.
salinity_perturbation = RandomPerturbations(-0.36696428571428574, 0.1)
set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
                                  salinity_perturbation)
## build the simulation
Δt = 1e-5
stop_time = 5 * 60
save_schedule = 5 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
