using Rasters, NCDatasets, CairoMakie
cab_file = "data/simulations/data/simulations/cabbeling_stepchange_dd_3min.nc"
σ = Raster(cab_file, name = :σ, lazy = true)
x, y, t = lookup(σ, :xC), lookup(σ, :zC), lookup(σ, :Ti)

yslice, snapshot = 62, 25
fig = Figure(size = (500, 500))
ax = Axis(fig[1, 1],
          xlabel = "x (m)",
          ylabel = "z (m)",
          title = "Potential density at t = $(round(t[snapshot] / 60)) minutes")
σ_snapshot = σ.data[:, yslice, :, snapshot]
hm = heatmap!(ax, x, z, σ_snapshot, colormap = :dense)
save("density.png", fig)
