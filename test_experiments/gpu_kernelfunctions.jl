using Printf
using Oceananigans: BuoyancyModels.ρ′, BuoyancyModels.get_temperature_and_salinity, BuoyancyModels.θ_and_sᴬ
using Oceananigans: BuoyancyModels.∂z_b
using Oceananigans: Operators.ℑzᵃᵃᶠ
import Oceananigans.BuoyancyModels.ρ′
Oceananigans.BuoyancyModels.ρ′(i, j, k, grid, eos, θ, sᴬ, pᵣ) = ρ′(θ_and_sᴬ(i, j, k, θ, sᴬ)..., pᵣ, eos)

## A test model on the GPU
architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))
DNS_resolution = (Nx = 10, Ny = 10, Nz = 100)

## Setup the model
@info "Model setup"
model = DNS(architecture, DOMAIN_EXTENT, DNS_resolution, diffusivities;
            reference_density = REFERENCE_DENSITY, zgrid_stretching = false)

## set initial conditions
@info "Setting initial conditions"
T₀ᵘ = -1.5
S₀ᵘ = 34.568
cabbeling = CabbelingUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
depth = find_depth(model, INTERFACE_LOCATION)
profile_function = StepChange(depth)

## Salinity noise
depths = find_depth(model, [INTERFACE_LOCATION + 0.02, INTERFACE_LOCATION - 0.02])
scales = similar(depths)
fill!(scales, 2e-4)
initial_noise = SalinityNoise(depths, scales)
@info "Building DNS"
dns = TwoLayerDNS(model, profile_function, initial_conditions; initial_noise)

@info "Setting two layer initial conditions"
set_two_layer_initial_conditions!(dns)

@inline function Kᵥ(i, j, k, grid, b::SeawaterBuoyancy, C, w)

    if ∂z_b(i, j, k, grid, b, C) != 0
        (-ℑzᵃᵃᶜ(i, j, k, w) * buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, b, C)) / ∂z_b(i, j, k, grid, b, C)
    else
        0
    end

end
function InferredVerticalDiffusivity(model)

    grid = model.grid
    b = model.buoyancy.model
    C = model.tracers
    w = model.velocities.w

    return KernelFunctionOperation{Center, Center, Center}(Kᵥ, grid, b, C, w)
end
