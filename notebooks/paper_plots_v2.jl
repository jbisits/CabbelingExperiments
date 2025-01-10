using GLMakie, JLD2, GibbsSeaWater
cd(@__DIR__)
dns_output = "../dns_runs/cabbeling_stepchange_nothing_780min/analysis_and_field_snapshots.jld2"
isothermal_output = "../dns_runs/isothermal_stepchange_nothing_780min/isothermal_analysis_and_field_snapshots.jld2"
## Figure theme
publication_theme = Theme(font="CMU Serif", fontsize = 20,
                          Axis=(titlesize = 22,
                                xlabelsize = 20, ylabelsize = 20,
                                xgridstyle = :dash, ygridstyle = :dash,
                                xtickalign = 0, ytickalign = 0,
                                yticksize = 6.5, xticksize = 6.5),
                          Legend=(framecolor = (:black, 0.5),
                                  backgroundcolor = (:white, 0.5),
                                  labelsize = 20),
                          Colorbar=(ticksize=12,
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

fig = Figure(size = (1400, 600), px_per_unit=16)

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

# Figure two, DNS states multiple panels
## Load data to create figure
file = jldopen(dns_output)

x = file["dims/x"]
y = file["dims/y"]
z = file["dims/z"]
times = file["dims/timestamps"]

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
z_xy_interface = [zmid, zmid + 3, 875, 925] # Need to save these midsections when gadi is up

snapshots = [0, 3, 13, 24]
slices = [(density_xz_face1 = file["σ"]["σ_xzslice/σ_$(60.0 * snapshot)"],
          density_middle_slice = file["σ"]["σ_xyslice/σ_$(60.0 * snapshot)"],
          velocity_yz = file["w"]["w_yzslice/w_$(60.0 * snapshot)"][:, 1:end-1],
          velocity_zmean = file["w"]["w_zmean/w_$(60.0 * snapshot)"][:, :, 1]) for snapshot ∈ snapshots]
close(file)

## One row
fig = Figure(size = (1600, 700), px_per_unit = 16)

ax = [Axis3(fig[1, i],
           aspect=(1/3, 1/3, 1),
           xlabel = "x (m)",
           ylabel = "y (m)",
           zlabel = "z (m)",
           xlabeloffset = 50,
           ylabeloffset = 50,
           zlabeloffset = 70,
           xlabelsize = 18,
           ylabelsize = 18,
           zlabelsize = 18,
           xticklabelsize = 18,
           yticklabelsize = 18,
           zticklabelsize = 18,
           zlabelrotation = π / 2,
           limits = ((x[1], x[end]), (y[1], y[end]), (z[1], z[end])),
           elevation = π / 6.5,
           azimuth = 1.25π,
           xspinesvisible = false,
           yspinesvisible = false,
           zspinesvisible = false,
           zgridvisible = false,
           protrusions = 20
           ) for i ∈ 1:4]

velocity_colorrange = maximum(abs.(extrema(slices[4].velocity_yz)))  * [-1, 1]
density_colorrange = extrema(slices[2].density_xz_face1)

for i ∈ eachindex(snapshots)
    sf_σ = surface!(ax[i], x_xz, y_xz_density, z_xz; color = slices[i].density_xz_face1,
                    colormap = :dense,
                    colorrange = density_colorrange
                    )
    sf_w = surface!(ax[i], x_xz_velocity, y_xz, z_yz; color = slices[i].velocity_yz, colormap = :balance,
                    colorrange = velocity_colorrange, backlight = 1f0, shading = FastShading)
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
rowgap!(fig.layout, 1, Relative(0.08))
fig
##
save("dns_schematic_ts_horizontal.png", fig)

## Extract dns ouput
file = jldopen(dns_output)
timestamps = file["dims/timestamps"]
zF = file["dims/zF"]
zC = file["dims/z"]

# Horizontally averaged profile
ha_σ = hcat([reshape(file["σ/σ_$(t)"], :) for t ∈ timestamps]...)
ρ_max = file["attrib/ρ_max"]

# Diffusivity
κₛ = reverse(file["diffusivity/κₛ"], dims = 1)
∫κₛ = file["diffusivity/∫κₛ"]

# Energetics
ρ₀ = file["attrib/ρ₀"]
g = file["attrib/g"]
 # the correction is due to incorrect scaling by area in the original computation
 bpe = file["energetics"]["∫Eb"]  * ρ₀ * (0.1^2 / 0.07^2)
  pe = file["energetics"]["∫Ep"]  * ρ₀
  ek = file["energetics"]["∫Ek"]  * ρ₀
   ε = file["energetics"]["∫ε"]   * ρ₀
∫gρw = file["energetics"]["∫gρw"] * ρ₀
  Δt = file["dims"]["Δt"]

ape = pe .- bpe
dₜek = diff(ek) ./ Δt
dₜpe = diff(pe) ./ Δt
dₜbpe = Φd = diff(bpe) ./ Δt
dₜape = diff(ape) ./ Δt
time_interp = 0.5 * (timestamps[1:end-1] .+ timestamps[2:end])

close(file)

# isothermal output
file = jldopen(isothermal_output)
κₛ_isothermal = file["diffusivity"]["κₛ"]
∫κₛ_isothermal = file["diffusivity"]["∫κₛ"]
close(file)
## Results figure:
#  HA profile evolution, HA effective diffusivity, HA volume integrated effective diffusivity
#  Energetics

fig = Figure(size = (1400, 1200), px_per_unit = 16)
time_interp_mins = time_interp ./ 60
time_ticks = (0:200:1080, string.(0:200:1080))
restricted_time = 1:600
# HA profile
ρ_anomaly = 0
snapshots = [1, 50, 100, 400, 1080] .+ 1
ax1 = Axis(fig[1, 1], xlabel = "σ₀ (kgm⁻³) ", ylabel = "z (m)",
          title = "(a) Horizontally averaged density profile")
lines!(ax1, ha_σ[:, 1] .- ρ_anomaly, zC, label = "Initial profile", color = :black)
vlines!(ax1, ρ_max - ρ_anomaly, color = :red, linestyle = :dash, label = "Maximum σ₀")
for (i, t) ∈ enumerate(snapshots)
    lines!(ax1, ha_σ[:, t] .- ρ_anomaly, zC, label = "$(timestamps[t] / 60) mins", alpha = 0.75)
end
axislegend(ax1, nbanks = 2, position = :lb, labelsize = 16)

# Diffusivity
ax2 = Axis(fig[1, 2], title = "(b) Horizontally averaged salinity effective diffusivity",
            xlabel = "time (mins)", ylabel = "z (m)", xticks = time_ticks)
hm = heatmap!(ax2, time_interp_mins[restricted_time], zC, log10.(abs.(κₛ[:, restricted_time]')),
                colorrange = (log10(1e-8), log10(1)), colormap = :tempo )
Colorbar(fig[1, 3], hm, label = "Effective diffusivity (m²s⁻¹, log10)")

ax4 = Axis(fig[2, 2], xlabel = "time (mins)", ylabel = "Effective diffusivity (m²s⁻¹, log10)",
            title = "(d) Depth integrated horizontally averaged\nsalinity effective diffusivity",
            xticks = time_ticks)
lines!(ax4, time_interp_mins[restricted_time], log10.(abs.(∫κₛ_isothermal[restricted_time])),
        label = "Isothermal", color = :cyan)
lines!(ax4, time_interp_mins[restricted_time], log10.(abs.(∫κₛ[restricted_time])),
        label = "Cabbeling", color = :pink)
hlines!(ax4, log10(1e-7), label = "Parameterised salinity diffusivity", linestyle = :dash,
        color = :black)
axislegend(ax4, position = :rt, labelsize = 17)

xlims!(ax4, 1e-5, restricted_time[end])
linkxaxes!(ax2, ax4)
hidexdecorations!(ax2, ticks = false)
hideydecorations!(ax2, ticks = false)
linkyaxes!(ax1, ax2)

# Energetics
ax3 = Axis(fig[2, 1], xlabel = "time (mins)", ylabel = "Watts",
          title = "(c) Time derivative of energetic quantities")
energy_scale = 1e8
lines!(ax3, time_interp_mins[restricted_time],  dₜpe[restricted_time] .* energy_scale,
        label = "dₜPE")
lines!(ax3, time_interp_mins[restricted_time], dₜbpe[restricted_time] .* energy_scale,
        label = "dₜBPE")
lines!(ax3, time_interp_mins[restricted_time], dₜape[restricted_time] .* energy_scale,
        label = "dₜAPE", color = :red)
# lines!(ax3, time_interp_mins[1:400],  dₜek[1:400] .* energy_scale, label = "dₜEk", color = :green)
Label(fig[2, 1, Top()], halign = :left, L"\times 10^{-8}", fontsize = 16)
axislegend(ax3, position = :rt, orientation = :horizontal)
# ylims!(ax3, -0.05, 0.4)

topspinecolor = leftspinecolor = rightspinecolor = bottomspinecolor = :gray
zoom_window = 1:100

ax_inset = Axis(fig[2, 1],
    width=Relative(0.68),
    height=Relative(0.572),
    halign=0.85,
    valign=0.6,
    backgroundcolor=:white,
    xticklabelsize = 12,
    yticklabelsize = 12;
    topspinecolor, leftspinecolor, rightspinecolor, bottomspinecolor
    )
ylims!(ax_inset, -0.05, 0.4)
lines!(ax_inset, time_interp_mins[zoom_window],  dₜpe[zoom_window] .* energy_scale)
lines!(ax_inset, time_interp_mins[zoom_window], dₜbpe[zoom_window] .* energy_scale)
lines!(ax_inset, time_interp_mins[zoom_window], dₜape[zoom_window] .* energy_scale, color = :red)
hlines!(ax_inset, 0, linestyle = :dash, color = :grey)
text!(ax3, 148, 3.85, text = L"\times 10^{-8}", fontsize = 12)
translate!(ax_inset.scene, 0, 0, 10)
# this needs separate translation as well, since it's drawn in the parent scene
translate!(ax_inset.elements[:background], 0, 0, 9)

using GLMakie.GeometryBasics
ϵ = 1e-12
rectangle_corners = Point2f[(0, minimum(dₜape)*energy_scale - ϵ),
                            (time_interp_mins[zoom_window[end]], minimum(dₜape)*energy_scale-ϵ),
                            (time_interp_mins[zoom_window[end]], maximum(dₜpe)*energy_scale+ϵ),
                            (0, maximum(dₜpe)*energy_scale + ϵ)]
poly!(ax3, rectangle_corners, strokecolor = :grey, strokewidth = 1, color = :transparent)
xs = [time_interp_mins[310]]
ys = [1]
us = [time_interp_mins[zoom_window[end]+10] - time_interp_mins[315]]
vs = [0]
arrows!(ax3, xs, ys, us, vs, color = :gray)
fig

##
save("results.png", fig)


## Energy budget
ε_interp = 0.5 * (ε[1:end-1] .+ ε[2:end])
∫gρw_interp = 0.5 * (∫gρw[1:end-1] .+ ∫gρw[2:end])
RHS = -ε_interp .- ∫gρw_interp
fig, ax = lines(time_interp, dₜek)
lines!(ax, time_interp, RHS)
fig

## energetic quantites
fig, ax = lines(timestamps, pe)
lines!(ax, timestamps, bpe)
fig

lines(dₜape[3:100])
findlast(dₜape[1:400] .> 0)
fig, ax = lines(time_interp, dₜpe)
lines!(ax, time_interp, dₜbpe)
fig

## Graphical abstract
fig = Figure(size = (1600, 700), px_per_unit = 16)

ax = [Axis3(fig[1, i],
           aspect=(1/3, 1/3, 1),
           xlabel = "x (m)",
           ylabel = "y (m)",
           zlabel = "z (m)",
           xlabeloffset = 50,
           ylabeloffset = 50,
           zlabeloffset = 70,
           xlabelsize = 18,
           ylabelsize = 18,
           zlabelsize = 18,
           xticklabelsize = 18,
           yticklabelsize = 18,
           zticklabelsize = 18,
           zlabelrotation = π / 2,
           limits = ((x[1], x[end]), (y[1], y[end]), (z[1], z[end])),
           elevation = π / 6.5,
           azimuth = 1.25π,
           xspinesvisible = false,
           yspinesvisible = false,
           zspinesvisible = false,
           zgridvisible = false,
           protrusions = 20
           ) for i ∈ 1:4]

velocity_colorrange = maximum(abs.(extrema(slices[4].velocity_yz)))  * [-1, 1]
density_colorrange = extrema(slices[2].density_xz_face1)

for i ∈ eachindex(snapshots)
    sf_σ = surface!(ax[i], x_xz, y_xz_density, z_xz; color = slices[i].density_xz_face1,
                    colormap = :dense,
                    colorrange = density_colorrange
                    )
    sf_w = surface!(ax[i], x_xz_velocity, y_xz, z_yz; color = slices[i].velocity_yz, colormap = :balance,
                    colorrange = velocity_colorrange, backlight = 1f0, shading = FastShading)
    surface!(ax[i], x, y, z_xy_top; color = slices[i].velocity_zmean, colormap = :balance,
            colorrange = velocity_colorrange
            )
    z_xy_mid = z[z_xy_interface[i]] * ones(length(x), length(y))
    surface!(ax[i], x, y, z_xy_mid; color = slices[i].density_middle_slice, colormap = :dense,
           colorrange = i == 1 ? extrema(slices[i].density_middle_slice) : density_colorrange,
            overdraw = true)
    hidedecorations!(ax[i])
end
# rowgap!(fig.layout, 1, Relative(0.08))
fig
##
save("graphical_abstract.png", fig)
