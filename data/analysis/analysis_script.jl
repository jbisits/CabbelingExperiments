using TwoLayerDirectNumericalShenanigans, NCDatasets, Rasters, CairoMakie

tracers = joinpath(@__DIR__, "tracers.nc")
computed_output = joinpath(@__DIR__, "computed_output.nc")

## Append output if it is not present
computed_variables = keys(NCDataset(computed_output))
if "∫κᵥ" ∉ computed_variables
    @info "Appeneding inferred vertical diffusivity"
    TLDNS.inferred_vertical_diffusivity!(computed_output, :∫ₐw′T′, :∫ₐ∂T∂z)
end
ds_attributes = keys(NCDataset(computed_output)ds.attrib)
if "λ_B" ∉ keys(ds_attributes)
    @info "Appending Kolmogorov and Batchelor scale"
    TLDNS.kolmogorov_and_batchelor_scale!(computed_output)
end

## Animate tracers
pred_Tₗ, pred_Sₗ = NCDataset(tracers) do ds
    (ds.attrib["Predicted equilibrium Tₗ"],
     ds.attrib["Predicted equilibrium Sₗ"])
end

@info "Reading into T into Raster"
T_rs = Raster(tracers, lazy = true, name = :T)
@info "Animating temperature"
colormap = cgrad(:thermal)[2:end-1]
colorrange = extrema(T_rs[:, :, :, 1])
lowclip = cgrad(:thermal)[1]
highclip = cgrad(:thermal)[end]
vline = pred_Tₗ
animate_2D_field(T_rs, 10, 10; colormap, colorrange, highclip, lowclip, vline)
@info "Reading into S into Raster"
S_rs = Raster(tracers, lazy = true, name = :S)
@info "Animating salinity"
colormap = cgrad(:haline)[2:end-1]
colorrange = extrema(S_rs[:, :, :, 1])
lowclip = cgrad(:haline)[1]
highclip = cgrad(:haline)[end]
vline = pred_Sₗ
animate_2D_field(S_rs, 10, 10; colormap, colorrange, highclip, lowclip, vline)

## Animate and plot computed output
pred_max_density = NCDataset(computed_output) do ds
    ds.attrib["Predicted maximum density"]
end
@info "Reading into σ into Raster"
σ_rs = Raster(computed_output, lazy = true, name = :σ)
@info "Animating density"
colormap = cgrad(:dense)[2:end-1]
colorrange = (minimum(σ_rs[:, :, :, 1]), pred_max_density)
lowclip = cgrad(:dense)[1]
highclip = cgrad(:dense)[end]
vline = pred_max_density
animate_2D_field(σ_rs , 10, 10; colormap, colorrange, highclip, lowclip, vline)

## Scalar diagnostics
plot_scalar_diagnostics(computed_output)

## Hovmoller of inferred vertical diffusivity
hovmoller(computed_output, "∫ₐκᵥ", unit = "m²s⁻¹")
