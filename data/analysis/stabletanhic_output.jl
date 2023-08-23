using DirectNumericalCabbelingShenanigans.OutputUtilities

## Initial snapshots
visualise_snapshot(T_ts, "Θ (°C)", 10, 10, 11)
visualise_snapshot(S_ts, "S (gkg⁻ꜝ)", 10, 10, 11; colormap = :haline)
visualise_snapshot(σ₀_ts, "σ₀ (kgm⁻³)", 10, 10, 10; colormap = :dense)

## Animations (x-z)
#  Temperatures
animate_2D_field(T_ts, "Θ (°C)", 10, 10)
## Salinity
animate_2D_field(T_ts, "S (gkg⁻¹)", 10, 10; colormap = :haline)
## Density (x-z)
animate_2D_field(σ₀_ts, "σ₀ (kgm⁻³)", 10, 10; colormap = :dense)
