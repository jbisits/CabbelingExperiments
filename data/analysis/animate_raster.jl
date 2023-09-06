using TwoLayerDirectNumericalShenanigans, CairoMakie, Rasters, NCDatasets

## Animations (x-z)
@info "Reading into T into Raster"
rs = Raster(stable_lp, name = :T)
@info "Animating temperature"
colormap = cgrad(:thermal)[2:end-1]
colorrange = extrema(rs)
lowclip = cgrad(:thermal)[1]
highclip = cgrad(:thermal)[end]
#  Temperatures
animate_2D_field(rs, 10, 10; colormap, colorrange, highclip, lowclip)
## Salinity
@info "Reading into S into Raster"
rs = Raster(stable_lp, name = :S)
@info "Animating salinity"
colormap = cgrad(:haline)[2:end-1]
colorrange = extrema(rs)
lowclip = cgrad(:haline)[1]
highclip = cgrad(:haline)[end]
animate_2D_field(rs, 10, 10; colormap, colorrange, highclip, lowclip)
## Density (x-z)
# @info "Reading into S into Raster"
# rs = Raster(stable_lp, name = :σ)
# @info "Animating density"
# colormap = cgrad(:dense)[2:end-1]
# colorrange = extrema(sr)
# lowclip = cgrad(:dense)[1]
# highclip = cgrad(:dense)[end]
# animate_2D_field(rs, 10, 10; colormap, colorrange, highclip, lowclip)
