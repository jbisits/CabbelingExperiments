using TwoLayerDirectNumericalShenanigans

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions
T₀ᵘ = -1.5
S = 16.3240507756
isohaline = IsohalineUpperLayerInitialConditions(S, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(isohaline)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 3500.0)

## `GaussianProfile` with `RandomPerturbations`
tracer_perturbation = SalinityGaussianProfile(INTERFACE_LOCATION, INTERFACE_LOCATION / 1.1,
                                        100.0, 10.0)
# z = znodes(model.grid, Center(), Center(), Center())
# depth_idx = findfirst(z .> INTERFACE_LOCATION / 1.1)
depth = find_depth(model, INTERFACE_LOCATION / 1.1)
initial_noise = SalinityNoise(-0.34077380952380953, 0.001)
dns = TwoLayerDNS(model, profile_function, initial_conditions;
                  tracer_perturbation, initial_noise)

set_two_layer_initial_conditions!(dns)


## build the simulation
Δt = 1e-5
stop_time = 10 * 60 # seconds (in simulation time)
save_schedule = 5 # seconds
simulation = DNS_simulation_setup(dns, Δt, stop_time, save_schedule)

## Run the simulation
run!(simulation)
