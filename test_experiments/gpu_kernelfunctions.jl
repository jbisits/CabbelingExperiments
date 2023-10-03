using Printf
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
"""
    DensityField(model, reference_pressure)
Return an `KernelFunctionOperation` at `(Center, Center, Center)` that computes the
potential density from the salinity and temperature tracers in `model` at `reference_pressure`.
"""
function DensityField(model)

    #tracers = (S = model.tracers.S, T = model.tracers.T)

   return KernelFunctionOperation{Center, Center, Center}(Oceananigans.BuoyancyModels.ρ′, model.grid, model.buoyancy.model.equation_of_state)

end
