## Script to plot a hovmoller of effective diffusivity and compute an average vertical
#  diffusivty.

using NCDatasets, CairoMakie, StatsBase

## Temperatre
# T_diff = "T_effective_diffusivity.nc"

# dsT = NCDataset(T_diff, "a")

# diffusivity = dsT[:κ_effectiveT][:, :]
# replace!(diffusivity, Inf => NaN)
# replace!(diffusivity, -Inf => NaN)
# t = dsT[:time_derivative][:]
# ez = dsT[:equivalent_z][:]

# fig = Figure(size = (500, 500))
# axT = Axis(fig[1, 1], xlabel = "time (s)", ylabel = "Equivalent z (m)",
#            title = "Effective temperature diffusivity")
# hm = heatmap!(axT, ez, t, log10.(diffusivity), colormap = :balance)
# Colorbar(fig[1, 2], hm, label = "Diffusivity (m²s⁻¹)")
# hidexdecorations!(axT, grid = false)
# @info "Saving figure"
# save("Tdiffusivity_hov.png", fig)

# mean_diff = "κ_effectiveT_mean"

# if mean_diff ∉ keys(dsT)

#     @info "Computing mean"
#     defVar(dsT, "κ_effectiveT_mean", Float64, tuple("time_derivative"),
#             attrib = Dict("long_name" => "Vertically averaged effective diffusivity"))

#     for (i, c) ∈ enumerate(eachcol(diffusivity))

#         dsT[mean_diff][i] = mean(c[.!isnan.(c)])

#     end

# end

# close(dsT)

## Salintiy
S_diff = "S_effective_diffusivity.nc"
dsS = NCDataset(S_diff, "a")

diffusivity = dsS[:κ_effectiveS][:, :]
replace!(diffusivity, Inf => NaN)
replace!(diffusivity, -Inf => NaN)
t = dsS[:time_derivative][:]
ez = dsS[:equivalent_z][:]

# fig = Figure(size = (500, 500))
# axS = Axis(fig[1, 1], xlabel = "time (s)", ylabel = "Equivalent z (m)",
#            title = "Effective Salinity diffusivity")

# hm = heatmap!(axS, ez, t, diffusivity, colormap = :balance)
# Colorbar(fig[1, 2], hm, label = "Diffusivity (m²s⁻¹)")

# save("Sdiffusivity_hov.png", fig)

mean_diff = "κ_effectiveS_mean"

if mean_diff ∉ keys(dsS)

    defVar(dsS, "κ_effectiveS_mean", Float64, tuple("time_derivative"),
            attrib = Dict("long_name" => "Vertically averaged effective diffusivity"))

    for (i, c) ∈ enumerate(eachcol(diffusivity))

        dsS[mean_diff][i] = mean(c[.!isnan.(c)])

    end

end

close(dsS)
