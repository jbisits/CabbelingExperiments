## Script to plot a hovmoller of effective diffusivity and compute an average vertical
#  diffusivty.

using NCDatasets, CairoMakie, StatsBase

## Temperatre
T_diff = "T_effective_diffusivity.nc"

fig = Figure(size = (800, 1000))
axT = Axis(fig[1, 1], xlabel = "time (s)", ylabel = "Equivalent z (m)",
           title = "Effective temperature diffusivity")

dsT = NCDataset(T_diff)

diffusivity = dsT[:κ_effectiveT][:, :]
replace!(diffusivity, Inf => NaN)
t = dsT[:time_derivative][:]
ez = dsT[:equivalent_z][:]
hm = heatmap!(axT, ez, t, diffusivity, colormap = :balance)
Colorbar(fig[1, 2], hm, colorscale = :log10, label = "Diffusivity (m²s⁻¹)")
hidexdecorations!(axT, grid = false)

mean_diff = "κ_effectiveT_mean"

if meand_diff ∉ keys(dsT)

    defVar(dsT, "κ_effectiveT_mean", Float64, tuple("time_derivative"),
            attrib = Dict("long_name" => "Vertically averaged effective diffusivity"))

    for (i, c) ∈ enumerate(eachcol(diffusivity))

        dsT[mean_diff][i] = mean(c[.!isnan.(c)])

    end

end

close(dsT)

## Salintiy
S_diff = "S_effective_diffusivity.nc"

axS = Axis(fig[2, 1], xlabel = "time (s)", ylabel = "Equivalent z (m)",
           title = "Effective Salinity diffusivity")

dsS = NCDataset(S_diff)

diffusivity = dsS[:κ_effectiveS][:, :]
replace!(diffusivity, Inf => NaN)
hm = heatmap!(axS, ez, t, diffusivity, colormap = :balance)
Colorbar(fig[2, 2], hm, colorscale = :log10, label = "Diffusivity (m²s⁻¹)")

save("diffusivity_hov.png", fig)

mean_diff = "κ_effectiveS_mean"

if meand_diff ∉ keys(dsS)

    defVar(dsS, "κ_effectiveS_mean", Float64, tuple("time_derivative"),
            attrib = Dict("long_name" => "Vertically averaged effective diffusivity"))

    for (i, c) ∈ enumerate(eachcol(diffusivity))

        dsS[mean_diff][i] = mean(c[.!isnan.(c)])

    end

end

close(dsS)
