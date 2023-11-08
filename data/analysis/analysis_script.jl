using TwoLayerDirectNumericalShenanigans, NCDatasets, CairoMakie

tracers = joinpath(@__DIR__, "tracers.nc")
computed_output = joinpath(@__DIR__, "computed_output.nc")

# Animate tracers
@info "Animating tracers"
TLDNS.animate_tracers(tracers)

# Animate density
@info "Animating density"
TLDNS.animate_density(computed_output)

# Scalar diagnostics
@info "Scalar diagnostics"
TLDNS.plot_scalar_diagnostics(computed_output)

# Hovmoller of inferred vertical diffusivity
@info "Hovmoller of horizontally integrated inferred vertical diffusivity"
TLDNS.hovmoller(computed_output, "∫ₐκᵥ", unit = "m²s⁻¹")
