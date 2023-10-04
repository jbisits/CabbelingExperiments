using Printf
using Oceananigans: BuoyancyModels.ρ′, BuoyancyModels.get_temperature_and_salinity, BuoyancyModels.θ_and_sᴬ
using Oceananigans: BuoyancyModels.∂z_b
using Oceananigans: Operators.ℑzᵃᵃᶠ
import Oceananigans.BuoyancyModels.ρ′
Oceananigans.BuoyancyModels.ρ′(i, j, k, grid, eos, θ, sᴬ, pᵣ) =
    ρ′(θ_and_sᴬ(i, j, k, θ, sᴬ)..., pᵣ, eos)

## Works on GPU.
# @inline function densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C)
#     T, S = get_temperature_and_salinity(b, C)
#     return  b.equation_of_state.reference_density + ρ′(i, j, k, grid, b.equation_of_state, T, S)
# end
# density(model) = density(model.buoyancy, model.grid, model.tracers)
# density(b, grid, tracers) = KernelFunctionOperation{Center, Center, Center}(densityᶜᶜᶜ, grid, b.model, tracers)
# DensityField(model) = Field(density(model))

## Computing at reference pressure (or reference geopotential height), not sure if this works on GPU
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
