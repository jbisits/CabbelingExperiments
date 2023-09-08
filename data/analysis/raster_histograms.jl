using RasterHistograms

TS_stack = RasterStack(sim_path, lazy = true, name = (:S, :T))
series_hist = slice(TS_stack, Ti())
stack_hist = RasterStackHistogram(TS_stack, (range(34.551, 34.7, step = 0.01), (range(-1.5, 0.5, step = 0.1))))
series_hist = RasterSeriesHistogram(series_hist, (range(34.551, 34.7, step = 0.01), (range(-1.5, 0.5, step = 0.1))))

fig = Figure(size = (500, 500))
ax = Axis(fig[1, 1];
          title = "Temperature and salinity joint distribution (unweighted)",
          xlabel = "Absolute salinity (gkg⁻¹)",
          ylabel = "Conservative temperature (°C)")
hm = heatmap!(ax, series_hist)
Colorbar(fig[1, 2], hm)
fig
