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

## Works on GPU.
# @inline function densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C)
#     T, S = get_temperature_and_salinity(b, C)
#     return  b.equation_of_state.reference_density + ρ′(i, j, k, grid, b.equation_of_state, T, S)
# end
# density(model) = density(model.buoyancy, model.grid, model.tracers)
# density(b, grid, tracers) = KernelFunctionOperation{Center, Center, Center}(densityᶜᶜᶜ, grid, b.model, tracers)
# DensityField(model) = Field(density(model))

## Computing at reference pressure (or reference geopotential height),
# works on GPU!!!!
@inline function densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C, parameters)
    T, S = get_temperature_and_salinity(b, C)
    pᵣ = parameters.pᵣ
    return  b.equation_of_state.reference_density + ρ′(i, j, k, grid, b.equation_of_state, T, S, pᵣ)
end
density(model, parameters) = density(model.buoyancy, model.grid, model.tracers, parameters)
density(b, grid, tracers, parameters) = KernelFunctionOperation{Center, Center, Center}(densityᶜᶜᶜ, grid, b.model, tracers, parameters)
DensityField(model, parameters) = Field(density(model, parameters))

parameters = (pᵣ = 0,)
d_field = DensityField(dns.model, parameters)
compute!(d_field)

## Center velocity `Field`
wᶜᶜᶜ(model) = KernelFunctionOperation{Center, Center, Center}(ℑzᵃᵃᶠ, model.grid, model.velocities.w)
w_center = Field(wᶜᶜᶜ(dns.model))
compute!(w_center)
## Center buoyancy perturbation field
b_field = BuoyancyField(model)
compute!(b_field)
## Face buoyancy gradient field
∂b∂z(model) = KernelFunctionOperation{Center, Center, Face}(∂z_b, model.grid, model.buoyancy, model.tracers)
b_vertical_grad = Field(∂b∂z(model))
compute!(b_vertical_grad)
Integral(b_field * w_center / b_vertical_grad)
