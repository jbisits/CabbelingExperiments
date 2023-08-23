using DirectNumericalCabbelingShenanigans.OutputUtilities

## Initial snapshots
visualise_snapshot(T_ts, "Θ (°C)", 10, 10, 11)
visualise_snapshot(S_ts, "S (gkg⁻ꜝ)", 10, 10, 11; colormap = :haline)
visualise_snapshot(σ₀_ts, "σ₀ (kgm⁻³)", 10, 10, 10; colormap = :dense)

## Animations (x-z)
#  Temperatures
animate_2D_field(T_ts, "Θ (°C)", 10, 10; colormap = cgrad(:thermal), colorrange = extrema(T_ts),
                                         highclip = cgrad(:thermal)[end], lowclip = cgrad(:thermal)[1])
## Salinity
animate_2D_field(T_ts, "S (gkg⁻¹)", 10, 10; colormap = cgrad(:haline), colorrange = extrema(S_ts),
                                            highclip = cgrad(:haline)[end], lowclip = cgrad(:haline)[1])
## Density (x-z)
animate_2D_field(σ₀_ts, "σ₀ (kgm⁻³)", 10, 10; colormap = cgrad(:dense), colorrange = extrema(σ₀_ts),
                                              highclip = cgrad(:dense)[end], lowclip = cgrad(:dense)[1])
