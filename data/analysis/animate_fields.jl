using CairoMakie
using DirectNumericalCabbelingShenanigans: animate_2D_field

## Animations (x-z)
@info "Animating temperature"
colormap = cgrad(:thermal)[2:end-1]
colorrange = extrema(T_ts[1])
lowclip = cgrad(:thermal)[1]
highclip = cgrad(:thermal)[end]
#  Temperatures
animate_2D_field(T_ts, "Θ (°C)", 10, 10; colormap, colorrange, highclip, lowclip)
## Salinity
@info "Animating salinity"
colormap = cgrad(:haline)[2:end-1]
colorrange = extrema(S_ts[1])
lowclip = cgrad(:haline)[1]
highclip = cgrad(:haline)[end]
animate_2D_field(S_ts, "S (gkg⁻¹)", 10, 10; colormap, colorrange, highclip, lowclip)
## Density (x-z)
@info "Animating density"
colormap = cgrad(:dense)[2:end-1]
colorrange = extrema(σ₀_ts[1])
lowclip = cgrad(:dense)[1]
highclip = cgrad(:dense)[end]
animate_2D_field(σ₀_ts, "σ₀ (kgm⁻³)", 10, 10; colormap, colorrange, highclip, lowclip)
