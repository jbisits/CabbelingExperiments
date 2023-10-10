using TwoLayerDirectNumericalShenanigans, CairoMakie, Rasters, NCDatasets

## Exract data from `.nc` file

pred_max_density = NCDataset(cab_noise) do ds
    ds["Predicted maximum density"]
end

## Animations (x-z)
@info "Reading into T into Raster"
T_rs = Raster(cab_noise, lazy = true, name = :T)
@info "Animating temperature"
colormap = cgrad(:thermal)[2:end-1]
colorrange = extrema(T_rs[:, :, :, 1])
lowclip = cgrad(:thermal)[1]
highclip = cgrad(:thermal)[end]
#  Temperatures
animate_2D_field(T_rs, 10, 10; colormap, colorrange, highclip, lowclip)
## Salinity
@info "Reading into S into Raster"
S_rs = Raster(cab_noise, lazy = true, name = :S)
@info "Animating salinity"
colormap = cgrad(:haline)[2:end-1]
colorrange = extrema(S_rs[:, :, :, 1])
lowclip = cgrad(:haline)[1]
highclip = cgrad(:haline)[end]
animate_2D_field(S_rs, 10, 10; colormap, colorrange, highclip, lowclip)
## Density (x-z)
@info "Reading into σ into Raster"
σ_rs = Raster(cab_noise, lazy = true, name = :σ)
@info "Animating density"
colormap = cgrad(:dense)[2:end-1]
colorrange = (minimum(σ_rs[:, :, :, 1]), pred_max_density)
lowclip = cgrad(:dense)[1]
highclip = cgrad(:dense)[end]
vline = pred_max_density
animate_2D_field(σ_rs , 10, 10; colormap, colorrange, highclip, lowclip)
