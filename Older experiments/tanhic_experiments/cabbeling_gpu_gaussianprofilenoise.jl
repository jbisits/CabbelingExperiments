using TwoLayerDirectNumericalShenanigans

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions
T₀ᵘ = -1.5
S₀ᵘ = 34.568
cabbeling = CabbelingUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 3500.0)

## `GaussianProfile` eith `RandomPerturbations`
tracer_perturbation = SalinityGaussianProfile(INTERFACE_LOCATION, INTERFACE_LOCATION / 1.1,
                                        100.0, 2.0)
# z = znodes(model.grid, Center(), Center(), Center())
# depth_idx = findfirst(z .> INTERFACE_LOCATION / 1.1)
depth = find_depth(model, INTERFACE_LOCATION / 1.1)
initial_noise = SalinityNoise(-0.34077380952380953, 0.1)
dns = TwoLayerDNS(model, profile_function, initial_conditions;
                  tracer_perturbation, initial_noise)

set_two_layer_initial_conditions!(dns)


## build the simulation
Δt = 1e-5
stop_time = 5 * 60
save_schedule = 5 # seconds
simulation = DNS_simulation_setup(dns, Δt, stop_time, save_schedule)

## Run the simulation
run!(simulation)
