using GLMakie, JLD2, GibbsSeaWater
cd(@__DIR__)

## Figure theme
publication_theme = Theme(font="CMU Serif", fontsize = 20,
                          Axis=(titlesize = 22,
                                xlabelsize = 20, ylabelsize = 20,
                                xgridstyle = :dash, ygridstyle = :dash,
                                xtickalign = 0, ytickalign = 0,
                                yticksize = 7.5, xticksize = 7.5),
                          Legend=(framecolor = (:black, 0.5),
                                  backgroundcolor = (:white, 0.5),
                                  labelsize = 20),
                          Colorbar=(ticksize=16,
                                    tickalign=1,
                                    spinewidth=0.5))
set_theme!(publication_theme)

## Figure one, schematic
S_star, Θ_star = 34.7, 0.5
S₀ᵘ = 34.58
Θ₀ᵘ = -1.5
slope = (Θ₀ᵘ - Θ_star) / (S₀ᵘ - S_star)
S_mix = range(S₀ᵘ, S_star, step = 0.000001)
Θ_mix = @. Θ₀ᵘ + (slope) * (S_mix - S₀ᵘ)
ρ_mix = gsw_rho.(S_mix, Θ_mix, 0)
ρ_max, ρ_max_idx = findmax(ρ_mix)
S_max, Θ_max = S_mix[ρ_max_idx], Θ_mix[ρ_max_idx]
Δρ_mix = ρ_max - gsw_rho(S_star, Θ_star, 0)

N = 1000
plotting_offset = 10
S_range, Θ_range = range(34.55, 34.705, length = N), range(-2, 1, length = N)
S_grid, Θ_grid = ones(N-(plotting_offset - 1)) .* S_range[plotting_offset:end]',
                 ones(N-(plotting_offset - 1))' .* Θ_range[plotting_offset:end]
ρ = gsw_rho.(S_grid, Θ_grid, 0)
ρ_star = gsw_rho(S_star, Θ_star, 0)
ρ_s = gsw_rho(S₀ᵘ, Θ₀ᵘ, 0)
find_Θ = findfirst(Θ_range .> -1.5)
find_S = findfirst(ρ[find_Θ, :] .> ρ_star)
S_iso, Θ_iso = S_range[find_S], Θ_range[find_Θ]
gsw_rho(S_iso, Θ_iso, 0)
αₗ, βₗ = gsw_alpha(S_star, Θ_star, 0), gsw_beta(S_star, Θ_star, 0)
m_initial = βₗ / αₗ
Θ_linear_initial = @. Θ_star + m_initial * (S_range - S_star)
αₘ, βₘ = gsw_alpha(S_max, Θ_max, 0), gsw_beta(S_max, Θ_max, 0)
m = βₘ / αₘ
Θ_linear = @. Θ_max + m * (S_range - S_max)

# Vertical profile
Nz = 200
z = range(-10, 0, length = 2*Nz)

Sˡ, Θˡ = fill(S_star, Nz), fill(Θ_star, Nz)
Sᵘ, Θᵘ = fill(S₀ᵘ, Nz), fill(Θ₀ᵘ, Nz)
S, Θ = vcat(Sˡ, Sᵘ), vcat(Θˡ, Θᵘ)
σ₀ = gsw_sigma0.(S, Θ)

Θ_mix_profile = (Θᵘ[1] + Θˡ[1]) / 2
S_mix_profile = (Sᵘ[1] + Sˡ[1]) / 2
mixing_interface = 190:210
S[mixing_interface] = fill(S_mix_profile, length(mixing_interface))
Θ[mixing_interface] = fill(Θ_mix_profile, length(mixing_interface))
σ₀_mix = gsw_sigma0.(S, Θ)

fig = Figure(size = (1000, 700))

ax = Axis(fig[1, 1],
          title = "(a) Two layer profile",
          xlabel = "σ₀",
          xaxisposition = :top,
          xlabelpadding = 2,
          ylabel = "z",
          yticksvisible = false,
          yticklabelsvisible = false,
          ygridvisible = false)
hidexdecorations!(ax, label = false)
hidespines!(ax)
lines!(ax, σ₀_mix, z, label = "Mixed water", color = :purple)
lines!(ax, σ₀[1:mixing_interface[1]-1], z[1:mixing_interface[1]-1], label = "Deep water", color = :red, linewidth = 2)
lines!(ax, σ₀[mixing_interface[end]+1:end], z[mixing_interface[end]+1:end], label = "Shallow water", color = :blue, linewidth = 2)
# density contours
vline_offset = 3
zvline = z[vline_offset-1:end-vline_offset]
Nvline = length(zvline)
lines!(ax, fill(σ₀_mix[end], Nvline), zvline, linestyle = :dot, linewidth = 2, color = (:blue, 0.5))
lines!(ax, fill(σ₀_mix[mixing_interface[10]], Nvline), zvline, linestyle = :dot, linewidth = 2, color = (:magenta, 0.5))
lines!(ax, fill(σ₀_mix[1], Nvline), zvline, linestyle = :dot, linewidth = 2, color = (:red, 0.5))
# drawn axis
arrows!(ax, [σ₀[end]-0.001], [z[end]], [1], [0], lengthscale = 0.011)
arrows!(ax, [σ₀[end]-0.001], [z[end]], [0], [-1], lengthscale = 10)
text!(ax, σ₀[end] + 0.0005, -2.5, text = "Cold/fresh\nwater",
      align = (:left, :center), color = :blue, justification = :center)
text!(ax, σ₀[1] + 0.0005, -7.5, text = "Warm/salty\nwater",
      align = (:left, :center), color = :red, justification = :center)
text!(ax, maximum(σ₀_mix)-0.0005, z[mixing_interface[11]], align = (:right, :center),
      text = "Mixed water", color = :purple)
text!(ax, σ₀_mix[end], z[10], align = (:left, :center),
      text = "σₒ shallow", color = :blue)
text!(ax, σ₀_mix[1], z[end-10], align = (:right, :center),
      text = "σₒ deep", color = :red)
text!(ax, maximum(σ₀_mix), z[end-10], align = (:right, :center),
      text = "σₒ maximum", color = :magenta)

ax2 = Axis(fig[1, 2];
          title = "(b) Cabbeling effect in S-Θ space",
          xlabel = "Absolute salinity",
          ylabel = "Conservative temperature",
          xticklabelsvisible = false,
          yticklabelsvisible = false,
          xgridvisible = false,
          xticksvisible = false,
          yticksvisible = false,
          ygridvisible = false,
          limits = (extrema(S_range), extrema(Θ_range)))
hidespines!(ax2)

# isopycnals and points
contour!(ax2, S_range[plotting_offset:end], Θ_range[plotting_offset:end], ρ'; levels = [ρ_s, ρ_star, ρ_max],
         color = [:blue, :red, :magenta], linestyle = :dot, linewidth = 2,
         labelsize = 18, label = "Isopycnals")
lines!(ax2, S_mix, Θ_mix, color = :purple, label = "Mixed water", linewidth = 0.8)
scatter!(ax2, [S_star], [Θ_star], color = :red, label = "Deep water")
scatter!(ax2, [S₀ᵘ], [Θ₀ᵘ], color = :blue, label = "Shallow water")
lines!(ax2, S_range[plotting_offset:end], Θ_linear_initial[plotting_offset:end],
        color = (:red, 0.5), label = "Tangent at\ndeep water", linestyle = :dashdot)
scatter!(ax2, S_max, Θ_max, color = :magenta, label = "Maximum\ndensity")

# wedge fill in
ρ_star_isopycnal = findall(ρ' .≈ ρ_star)
S_idx = [ρ_star_isopycnal[i][1] for i ∈ eachindex(ρ_star_isopycnal)]
Θ_idx = [ρ_star_isopycnal[i][2] for i ∈ eachindex(ρ_star_isopycnal)]
band_1_range = plotting_offset:S_idx[1]
S_band_1 = S_range[band_1_range]
Θ_lower_band_1 = fill(Θ_range[plotting_offset-1], length(band_1_range))
Θ_upper_band_1 = Θ_linear_initial[band_1_range]
band!(ax2, S_band_1, Θ_lower_band_1, Θ_upper_band_1, color = :red, alpha = 0.1)

band_2_range = S_idx
S_band_2 = S_range[band_2_range]
Θ_lower_band_2 = Θ_range[Θ_idx]
Θ_upper_band_2 = Θ_linear_initial[band_2_range]
band!(ax2, S_band_2, Θ_lower_band_2, Θ_upper_band_2, color = :red, alpha = 0.1)

# axis and text
arrows!(ax2, [S_range[10]], [Θ_range[5]], [1], [0], lengthscale = 0.15)
arrows!(ax2, [S_range[10]], [Θ_range[5]], [0], [1], lengthscale = 2.8)
arrows!(ax2, [S_mix[40000]], [Θ_mix[20000]], [0], [1], lengthscale = 0.30, color = :purple)
arrows!(ax2, [34.605], [-1.5], [-1], [0], lengthscale = 0.018, color = :red)
# Legend(fig[2, 2], ax2, orientation = :horizontal, nbanks = 3, tellheight = false)
text!(ax2, S_star-0.002, Θ_star, align = (:right, :bottom), text = "Deep water", color = :red)
text!(ax2, S₀ᵘ-0.001, Θ₀ᵘ, align = (:right, :bottom),text = "Shallow water", color = :blue)
text!(ax2, S_mix[40000], Θ_mix[20000], align = (:center, :top), text = "Mixed water", color = :purple)
text!(ax2, S_max, Θ_max, align = (:left, :top), text = "Maximum density", color = :magenta)
text!(ax2, 34.583, -1.7, align = (:left, :top),text = "σ₀ maximum", color = :magenta)
text!(ax2, 34.568, -1.7, align = (:right, :top),text = "σ₀ shallow", color = :blue)
text!(ax2, 34.605, -1.45, align = (:left, :top), text = "σ₀ deep", color = :red)
text!(ax2, 34.585, -1, align = (:right, :bottom),text = "Tangent to isopycnal\nat deep water", color = :red)

fig
colsize!(fig.layout, 1, Relative(0.3))
fig

##
save("schematic_2panel.png", fig)

## Figure two, DNS states multiple panels
## Load data to create figure
file_path = joinpath(@__DIR__, "density_and_w_snapshots.jld2")
file = jldopen(file_path)

x = file["σ/x"]
y = file["σ/y"]
z = file["σ/z"]
times = file["σ/times"]

x_xz = repeat(x, 1, length(z))
x_xy = repeat(x, 1, length(y))
y_xz = repeat(y, 1, length(z))
y_xy = repeat(y, 1, length(x))
z_xz = repeat(reshape(z, 1, length(z)), length(x), 1)
z_yz = repeat(reshape(z, 1, length(z)), length(x), 1)

x_xz_velocity = x[1] * ones(length(x), length(z))
x_xz_velocity_end = x[end] * ones(length(x), length(z))
y_xz_density = y[1] * ones(length(x), length(z))
y_xz_density_end = y[end] * ones(length(x), length(z))

z_xy_top = z[end] * ones(length(x), length(y))
# z_xy_mid = z[Int(length(z)/2)] * ones(length(x), length(y))
zmid = Int(length(z)/2)
z_xy_interface = [zmid, zmid + 3, 755, 780] # Need to save these midsections when gadi is up

snapshots = [0, 3, 12, 18]
slices = [(density_xz_face1 = file["σ/σ_xzslice_face1/σ_$(60.0 * snapshot)"],
          density_xz_faceend = file["σ/σ_xzslice_faceend/σ_$(60.0 * snapshot)"],
          density_xz_mid = file["σ/σ_xzslice/σ_$(60.0 * snapshot)"],
          density_middle_slice = snapshot < 12 ? file["σ/σ_xyslice/σ_$(60.0 * (snapshot))"] :
                                                 file["σ/σ_xyslice_migrate/σ_$(60.0 * (snapshot))"],
        #   density_middle_slice = file["σ/σ_xyslice/σ_$(60.0 * snapshot)"],
          velocity_yz = file["w/w_yzslice/w_$(60.0 * snapshot)"][:, 1:end-1],
          velocity_yz_face1 = file["w/w_yzslice_face1/w_$(60.0 * snapshot)"][:, 1:end-1],
          velocity_yz_faceend = file["w/w_yzslice_faceend/w_$(60.0 * snapshot)"][:, 1:end-1],
          velocity_zmean = file["w/w_zmean/w_$(60.0 * snapshot)"][:, :, 1]) for snapshot ∈ snapshots]
close(file)

## One row
fig = Figure(size = (1600, 800), fontsize = 18)

ax = [Axis3(fig[1, i],
           aspect=(1/3, 1/3, 1),
           xlabel = "x (m)",
           ylabel = "y (m)",
           zlabel = "z (m)",
           xlabeloffset = 50,
           ylabeloffset = 50,
           zlabeloffset = 70,
        #    xlabelsize = 14,
        #    ylabelsize = 14,
        #    zlabelsize = 14,
        #    xticklabelsize = 13,
        #    yticklabelsize = 13,
        #    zticklabelsize = 13,
           zlabelrotation = π / 2,
           limits = ((x[1], x[end]), (y[1], y[end]), (z[1], z[end])),
           elevation = π / 6.5,
           azimuth = 1.25π,
           xspinesvisible = false,
           yspinesvisible = false,
           zspinesvisible = false,
           zgridvisible = false,
           protrusions = 20,
           # persepectiveness = 0.1
           ) for i ∈ 1:4]

# velocity_colorrange = maximum([maximum(abs.(extrema(slices[i].velocity_yz_face1))) for i ∈ 1:4])  * [-1, 1]
velocity_colorrange = maximum(abs.(extrema(slices[4].velocity_yz_face1)))  * [-1, 1]
# velocity_colorrange = maximum(abs.(extrema(slices[4].velocity_zmean)))  * [-1, 1]
# density_colorrange = [minimum([minimum(slices[i].density_xz_face1) for i ∈ 1:4]),
#                       maximum([maximum(slices[i].density_xz_face1) for i ∈ 1:4])]
density_colorrange = extrema(slices[2].density_xz_face1)

for i ∈ eachindex(snapshots)
    sf_σ = surface!(ax[i], x_xz, y_xz_density, z_xz; color = slices[i].density_xz_face1,
                    colormap = :dense,
                    colorrange = density_colorrange
                    )
    sf_w = surface!(ax[i], x_xz_velocity, y_xz, z_yz; color = slices[i].velocity_yz_face1, colormap = :balance,
                    colorrange = velocity_colorrange, backlight = 1f0, shading = FastShading)
    sf_σ = surface!(ax[i], x_xz, y_xz_density_end, z_xz; color = fill(minimum(slices[i].density_xz_face1), size(slices[i].density_xz_faceend)),
                    # color = slices[i].density_xz_faceend,
                    colormap = :dense,
                    colorrange = density_colorrange,
                    )
    sf_w = surface!(ax[i], x_xz_velocity_end, y_xz, z_yz;
                    # color = slices[i].velocity_yz_faceend,
                    colormap = :balance,
                    colorrange = velocity_colorrange,
                    ambient = (0.85, 0.85, 0.85), backlight = 1f0)
    surface!(ax[i], x, y, z_xy_top; color = slices[i].velocity_zmean, colormap = :balance,
            colorrange = velocity_colorrange
            )
    z_xy_mid = z[z_xy_interface[i]] * ones(length(x), length(y))
    surface!(ax[i], x, y, z_xy_mid; color = slices[i].density_middle_slice, colormap = :dense,
           colorrange = i == 1 ? extrema(slices[i].density_middle_slice) : density_colorrange,
            overdraw = true)
    ax[i].title = "t = $(times[snapshots[i] + 1] / 60) minutes"
    if i == 1
        Colorbar(fig[2, 1:2], sf_σ, label = "σ₀ (kgm⁻³)", vertical = false, flipaxis = false)
        Colorbar(fig[2, 3:4], sf_w, label = "w (ms⁻¹)", vertical = false, flipaxis = false)
    end
    if i > 1
        hidezdecorations!(ax[i], ticks = false)
    end
end
# colgap!(fig.layout, 1, Relative(0.1))
# colgap!(fig.layout, 2, Relative(0.1))
# colgap!(fig.layout, 3, Relative(0.1))
rowgap!(fig.layout, 1, Relative(0.04))
fig
##
save("dns_schematic_ts_horizontal.png", fig)

## Horizontally averaged hovmollers, vertical velocity look suspicious
file = jldopen("ha_profile_snapshots.jld2")
timestamps = file["time"]
zF = file["zF"]
zC = file["zC"]
ρ_max = file["ρ_predicted_max"]

ha_σ = hcat([reshape(file["σ/σ_$(t)"], :) for t ∈ timestamps]...)
ha_w = hcat([reshape(file["w/w_$(t)"], :) for t ∈ timestamps]...)
close(file)

# iso = load("isothermal_fluxes_and_diff_longer_run.jld2")
iso = load("isothermal_fluxes_and_diff.jld2")
Δz_iso = diff(iso["z"])
replace!(iso["κₛ"], Inf => NaN)
replace!(iso["κₛ"], 0 => NaN)
replace!(iso["κₛ"], -Inf => NaN)
reverse!(iso["κₛ"], dims = 1)
κ_iso = similar(iso["∫κₛ"])
i = 1
for c ∈ eachcol(iso["κₛ"])
    find = findall(.!isnan.(c))
    κ_iso[i] = sum(c[find] .* Δz_iso[find]) / sum(Δz_iso[find])
    i += 1
end

cab = load("cabbeling_fluxes_and_diff_longer_run.jld2")
# cab = load("cabbeling_fluxes_and_diff.jld2")
cab_time = cab["time"]
Δz_cab = diff(cab["z"])
replace!(cab["κₛ"], Inf => NaN)
replace!(cab["κₛ"], 0 => NaN)
replace!(cab["κₛ"], -Inf => NaN)
reverse!(cab["κₛ"], dims = 1)

κ_cab = similar(cab["∫κₛ"])
j = 1
for c ∈ eachcol(cab["κₛ"])
    find = findall(.!isnan.(c))
    κ_cab[j] = sum(c[find] .* Δz_cab[find]) / sum(Δz_cab[find])
    j += 1
end

## energetics
energetics = load("cabbeling_energetics.jld2")
ρ₀ = energetics["ρ₀"]
pe = energetics["∫Ep_zref0"]*ρ₀
pe_anomaly = energetics["∫Ep_zref0_density_anomaly"]*ρ₀
bpe = energetics["∫Ebz✶_zref0"]*ρ₀
bpe_anomaly = energetics["∫Eb′z✶"]*ρ₀
ε = energetics["∫ε"]*ρ₀
ek = energetics["∫Ek"]*ρ₀
g = energetics["g"]
z = energetics["z"]
times = energetics["time"]

ape = pe .- bpe
ape_anomaly = pe_anomaly .- bpe_anomaly
dₜek = diff(ek) ./ diff(times)
dₜpe = diff(pe) ./ diff(times)
dₜbpe = Φd = diff(bpe) ./ diff(times)
dₜape = diff(ape) ./ diff(times)
dₜpe_anomaly = diff(pe_anomaly) ./ diff(times)
dₜbpe_anomaly = Φd = diff(bpe_anomaly) ./ diff(times)
dₜape_anomaly = diff(ape_anomaly) ./ diff(times)
time_interp = 0.5 * (times[1:end-1] .+ times[2:end])

## Results figure:
#  HA profile evolution, HA effective diffusivity, HA volume integrated effective diffusivity
#  Energetics

fig = Figure(size = (1200, 1200))
time_interp_mins = time_interp ./ 60
time_ticks = (0:200:1200, string.(0:200:1200))
# HA profile
ρ_anomaly = 0
snapshots = [1, 50, 100, 473, 1200] .+ 1
ax1 = Axis(fig[1, 1], xlabel = "σ₀ (kgm⁻³) ", ylabel = "z (m)",
          title = "(a) Horizontally averaged density profile", subtitle = "Cabbeling")
lines!(ax1, ha_σ[:, 1] .- ρ_anomaly, zC, label = "Initial profile", color = :black)
vlines!(ax1, ρ_max - ρ_anomaly, color = :red, linestyle = :dash, label = "Predicted maximum σ₀")
for (i, t) ∈ enumerate(snapshots)
    lines!(ax1, ha_σ[:, t] .- ρ_anomaly, zC, label = "$(timestamps[t] / 60) mins", alpha = 0.75)
end
axislegend(ax1, nbanks = 2, position = :lb, labelsize = 12)

# Diffusivity
ax2 = Axis(fig[1, 2], title = "(b) Horizontally averaged effective diffusivity", subtitle = "Cabbeling",
            xlabel = "time (mins)", ylabel = "z (m)", xticks = time_ticks)
hm = heatmap!(ax2, time_interp_mins, cab["z"], log10.(abs.(cab["κₛ"]')),
                colorrange = (log10(1e-8), log10(1)), colormap = :tempo )
Colorbar(fig[1, 3], hm, label = "Effective diffusivity (m²s⁻¹, log10)")

ax4 = Axis(fig[2, 2], xlabel = "time (mins)", ylabel = "Effective diffusivity (m²s⁻¹, log10)",
            title = "(d) Depth integrated horizontally averaged effective diffusivity",
            xticks = time_ticks)
lines!(ax4, time_interp_mins[1:660], log10.(abs.(κ_iso)), label = "Isothermal", color = :cyan)
lines!(ax4, time_interp_mins, log10.(abs.(κ_cab)), label = "Cabbeling", color = :pink)
hlines!(ax4, log10(1e-7), label = "Parameterised salinity diffusivity", linestyle = :dash, color = :black)
axislegend(ax4, position = :rt)

xlims!(ax4, 1e-5, 1200)
linkxaxes!(ax2, ax4)
hidexdecorations!(ax2, ticks = false)
hideydecorations!(ax2, ticks = false)

linkyaxes!(ax1, ax2)

# Energetics
ax3 = Axis(fig[2, 1], xlabel = "time (mins)", ylabel = "Watts",
          title = "(c) Time derivative of energetic quantities",
          subtitle = "Cabbeling",
          xticks = time_ticks)
energy_scale = 1e8
lines!(ax3, time_interp_mins, dₜpe_anomaly  .* energy_scale, label = "dₜPE")
lines!(ax3, time_interp_mins, dₜbpe_anomaly .* energy_scale, label = "dₜBPE")
lines!(ax3, time_interp_mins, dₜape_anomaly .* energy_scale, label = "dₜAPE", color = :red)
lines!(ax3, time_interp_mins, dₜek          .* energy_scale, label = "dₜEk", color = :green)
Label(fig[2, 1, Top()], halign = :left, L"\times 10^{-8}")
axislegend(ax3, position = :rt, orientation = :horizontal)

topspinecolor = leftspinecolor = rightspinecolor = bottomspinecolor = :gray
# bbox = BBox(275, 730, 230, 520)
# ax_inset = Axis(fig; bbox,  xticklabelsize = 12,
#             yticklabelsize = 12, topspinecolor, leftspinecolor, rightspinecolor, bottomspinecolor)
zoom_window = 1:220

ax_inset = Axis(fig[2, 1],
    width=Relative(0.67),
    height=Relative(0.5),
    halign=0.85,
    valign=0.6,
    backgroundcolor=:white,
    xticklabelsize = 12,
    yticklabelsize = 12;
    topspinecolor, leftspinecolor, rightspinecolor, bottomspinecolor
    )
lines!(ax_inset, time_interp_mins[zoom_window], dₜpe_anomaly[zoom_window]  .* energy_scale)
lines!(ax_inset, time_interp_mins[zoom_window], dₜbpe_anomaly[zoom_window] .* energy_scale)
lines!(ax_inset, time_interp_mins[zoom_window], dₜape_anomaly[zoom_window] .* energy_scale, color = :red)
lines!(ax_inset, time_interp_mins[zoom_window], dₜek[zoom_window]          .* energy_scale, color = :green)
hlines!(ax_inset, 0, linestyle = :dash)
text!(ax3, 300, 5.5, text = L"\times 10^{-8}", fontsize = 12)
translate!(ax_inset.scene, 0, 0, 10)
# this needs separate translation as well, since it's drawn in the parent scene
translate!(ax_inset.elements[:background], 0, 0, 9)

using GLMakie.GeometryBasics
ϵ = 1e-12
rectangle_corners = Point2f[(0, minimum(dₜape_anomaly)*energy_scale - ϵ),
                            (time_interp_mins[zoom_window[end]], minimum(dₜape_anomaly)*energy_scale-ϵ),
                            (time_interp_mins[zoom_window[end]], maximum(dₜpe_anomaly)*energy_scale+ϵ),
                            (0, maximum(dₜpe_anomaly)*energy_scale + ϵ)]
poly!(ax3, rectangle_corners, strokecolor = :grey, strokewidth = 1, color = :transparent)
xs = [time_interp_mins[310]]
ys = [1.1]
us = [time_interp_mins[zoom_window[end]+10] - time_interp_mins[300]]
vs = [0]
arrows!(ax3, xs, ys, us, vs, color = :gray)
fig

##
save("results.png", fig)

## Optimally mixed profile

## Horizontally averaged hovmollers, vertical velocity look suspicious
using Printf
file = jldopen("ha_profile_snapshots.jld2")
timestamps = file["time"]
zF = file["zF"]
zC = file["zC"]
ρ_max = file["ρ_predicted_max"]

ha_σ = hcat([reshape(file["σ/σ_$(t)"], :) for t ∈ timestamps]...)
close(file)

n = Observable(1)
σ_ha_plot = @lift ha_σ[:, $n] .- 1000
σ_ha_sorted_plot = @lift sort(ha_σ[:, $n], rev = true) .- 1000
time_title = @lift @sprintf("t=%1.2f minutes", timestamps[$n] / 60)

##
fig = Figure(size = (500, 500))
ax = Axis(fig[1, 1], title = time_title, xlabel = "σ₀ (kgm⁻³)", ylabel = "z (m)")
lines!(ax, ha_σ[:, 1].-1000, zC, label = "Initial density profile", color = :grey)
fig

N = length(zC)
σ_optimally_mixed = vcat(fill(ρ_max, N))
not_mixed = round(Int, N/2 - upper_proportion * N/2) # `upper_proportion` from amount required calculated below
σ_optimally_mixed[N-not_mixed:N] .= ha_σ[end, 1]
lines!(ax, σ_optimally_mixed, zC, label = "Optimally mixed", color = :black)

fig

ẑ = upper_proportion / lower_proportion - 1
hlines!(ax, ẑ, color = :orange)
## animate the ha profile onto the axis

lines!(ax, σ_ha_plot, zC, label = "In situ profile")
lines!(ax, σ_ha_sorted_plot, zC, label = "Sorted profile")

axislegend(ax, position = :lb)

frames = eachindex(timestamps)
record(fig, joinpath(pwd(), "ha_profiles.mp4"), frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")
    n[] = i
end

## Hovmoller
fig, ax, hm = heatmap(timestamps ./ 60, zC, ha_σ', colormap = :dense)
Colorbar(fig[1, 2], hm)
fig

ẑ = (upper_proportion * 0.5 + 0.5) - 1
ẑ = upper_proportion / lower_proportion - 1
hlines!(ax, ẑ, color = :orange)
fig

## Compare the PE values, proportional to g
ẑ = reverse(abs.(zC))
fig, ax = lines(ha_σ[:, 1], ẑ, label = "Initial profile")
lines!(ax, σ_optimally_mixed, ẑ, label = "Optimally mixed")
axislegend(ax)
fig
Ep = sum(ha_σ[:, 1] .* ẑ) * Δz
Ep_cabbeling_stable = sum(σ_optimally_mixed .* ẑ) * Δz
Ep < Ep_cabbeling_stable

## Optimally mixed
using SeawaterPolynomials.TEOS10
using SeawaterPolynomials: ρ
M = 5000
Sᵤ, Sₗ = 34.58, 34.7
Θᵤ, Θₗ = -1.5, 0.5
S_range = range(Sᵤ, Sₗ, length = M)
slope = (Θₗ - Θᵤ) / (Sₗ - Sᵤ)
Θ_linear =  @. Θₗ + slope * (S_range - Sₗ)

eos_vec = fill(TEOS10EquationOfState(reference_density = ρ₀), length(Θ_linear))
ρ_mix = ρ.(Θ_linear, S_range, 0, eos_vec)
ρmax, idx_ρmax = findmax(ρ_mix)
Smax, Θmax = S_range[idx_ρmax], Θ_linear[idx_ρmax]

initial_S_T = [Sᵤ Sₗ; Θᵤ Θₗ]
upper_proportion, lower_proportion = initial_S_T \ [Smax, Θmax]

N = 700
Sₗmix = fill(Smax, N)
Θₗmix = fill(Θmax, N)
upper_amount_mixed = round(Int, N * upper_proportion)
Sᵤmix = fill(Smax, upper_amount_mixed)
Θᵤmix = fill(Θmax, upper_amount_mixed)
Sᵤ_reduced = fill(34.58, N-upper_amount_mixed)
Θᵤ_reduced = fill(-1.5, N-upper_amount_mixed)

S_optimally_mixed = vcat(Sₗmix, Sᵤmix, Sᵤ_reduced)
Θ_optimally_mixed = vcat(Θₗmix, Θᵤmix, Θᵤ_reduced)

ρ_optimally_mixed = gsw_rho.(S_optimally_mixed, Θ_optimally_mixed, 0)

lines!(ax, ρ_optimally_mixed, zC)

## Where does dₜape remain less than zero
start_time = 472
fig, ax = lines(time_interp[start_time:end], dₜape_anomaly[start_time:end])
@info "After t = $(time_interp[start_time+1] ./60) dₜape < 0."

# Better way
dt_sign_idx = findfirst(reverse(dₜape_anomaly) .> 0)
reverse(time_interp)[dt_sign_idx-1] ./60

## Compare optimally mixed profile with profile when APE is negative definite
N = length(zC)
σ_optimally_mixed = vcat(fill(ρ_max - 1020, N))
not_mixed = round(Int, 0.57 * N/2)
σ_optimally_mixed[N-not_mixed:N] .= ha_σ[end, 1] .- 1020

fig, ax = lines(σ_optimally_mixed, zC)
lines!(ax, ha_σ[:, 473] .- 1020, zC)
fig

## And the same profile but Sorted, somthing like this could be useful
N = length(zC)
σ_optimally_mixed = vcat(fill(ρ_max, N))
not_mixed = round(Int, 0.57 * N/2)
σ_optimally_mixed[N-not_mixed:N] .= ha_σ[end, 1]

fig, ax = lines(σ_optimally_mixed, zC)
lines!(ax, sort(ha_σ[:, 473], rev = true), zC)
fig


## lengths for resolution
# η                    = 0.0019484207282447232
# λ_b                  = 0.0006161447341537292
# Δz                   = 0.0007142857142857784
# Δx                   = 0.0008064516129032279
