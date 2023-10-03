using Printf
using Oceananigans: BuoyancyModels.ρ′, BuoyancyModels.get_temperature_and_salinity
# import Oceananigans.BuoyancyModels.ρ′
# @inline ρ′(i, j, k, grid, eos, θ, sᴬ) = ρ′(θ_and_sᴬ(i, j, k, θ, sᴬ)..., 0, eos)
"""
    C(i, j, k, grid, C)
Get tracer `C` values for use in other function. There may be another way to do this for
`KernelFunctionOperation`s but I have not found it so will use this for now. **Note:** this
return the value of the tracer with no interpolation so if the tracer `C` is at
`(Center, Center, Center)` the value extracted will be at (Center, Center, Center)`.
"""
C(i, j, k, grid, C) = C[i, j, k]
"""
    σ(i, j, k, grid, model, reference_pressure)
Compute potential density `σ` at `reference_pressure` from salinity and temperature tracers
in `model`.
"""
σ(i, j, k, grid, tracers) = gsw_rho(C(i, j, k, grid, tracers.S),
                                    C(i, j, k, grid, tracers.T),
                                    0)

@inline function densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C)
    T, S = get_temperature_and_salinity(b, C)
    return  b.equation_of_state.reference_density + ρ′(i, j, k, grid, b.equation_of_state, T, S)
end
density(model) = density(model.buoyancy, model.grid, model.tracers)
density(b, grid, tracers) = KernelFunctionOperation{Center, Center, Center}(densityᶜᶜᶜ, grid, b.model, tracers)
DensityField(model) = Field(density(model))
"""
    DensityField(model, reference_pressure)
Return an `KernelFunctionOperation` at `(Center, Center, Center)` that computes the
potential density from the salinity and temperature tracers in `model` at `reference_pressure`.
"""
# function DensityField(model)

#     #tracers = (S = model.tracers.S, T = model.tracers.T)
#     parameters = (eos = model.buoyancy.model.equation_of_state)
#     return KernelFunctionOperation{Center, Center, Center}(densityᶜᶜᶜ, model.grid, parameters)

# end
