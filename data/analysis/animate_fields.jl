using DirectNumericalCabbelingShenanigans.OutputUtilities

## Animations (x-z)
@info "Animating temperature"
colormap = :thermal
colorrange = extrema(T_ts)
lowclip = :black
highclip = :yellow
#  Temperatures
animate_2D_field(T_ts, "Θ (°C)", 10, 10; colormap, colorrange)
## Salinity
@info "Animating salinity"
colormap = :haline
colorrange = extrema(S_ts)
lowclip = :black
highclip = :yellow
animate_2D_field(S_ts, "S (gkg⁻¹)", 10, 10; colormap, colorrange,
                                            highclip, lowclip)
## Density (x-z)
@info "Animating density"
colormap = cgrad(:dense)[2:end-1]
colorrange = extrema(σ₀_ts),
highclip = cgrad(:dense)[end]
lowclip = cgrad(:dense)[1]
animate_2D_field(σ₀_ts, "σ₀ (kgm⁻³)", 10, 10; colormap, colorrange,
                                              highclip, lowclip)
