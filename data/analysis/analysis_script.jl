using TwoLayerDirectNumericalShenanigans, NCDatasets, CairoMakie

tracers = joinpath(@__DIR__, "tracers.nc")
computed_output = joinpath(@__DIR__, "computed_output.nc")

# Tracer animations
@info "Animating tracers"
TLDNS.animate_tracers(tracers)

@info "Animating tracer distributions"
TLDNS.animate_tracer_distributions(tracers)

# Density animations
@info "Animating density"
TLDNS.animate_density(computed_output, "σ")

@info "Animating density distributions"
TLDNS.animate_density_distribution(computed_output)

# Scalar diagnostics
@info "Scalar diagnostics"
TLDNS.plot_scalar_diagnostics(computed_output)

# Hovmoller of inferred vertical diffusivity
@info "Hovmoller of horizontally integrated inferred vertical diffusivity"
TLDNS.hovmoller(computed_output, "∫ₐκᵥ", unit = "m²s⁻¹")
