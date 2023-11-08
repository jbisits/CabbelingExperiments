using TwoLayerDirectNumericalShenanigans, NCDatasets, CairoMakie

tracers = joinpath(@__DIR__, "tracers.nc")
computed_output = joinpath(@__DIR__, "computed_output.nc")

# Animate tracers
@info "Animating tracers"
animate_tracers(tracers)

# Animate density
@info "Animating density"
animate_density(computed_output)

# Scalar diagnostics
@info "Scalar diagnostics"
plot_scalar_diagnostics(computed_output)

# Hovmoller of inferred vertical diffusivity
@info "Hovmoller of horizontally integrated inferred vertical diffusivity"
hovmoller(computed_output, "∫ₐκᵥ", unit = "m²s⁻¹")
