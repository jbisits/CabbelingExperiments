# This version is the plots that have been updated during first revisions of the MS.
# Figures 1 and 2 are unchanged (except for labels). Figure 3, the previous results figure,
# will be broken up (likely into density evolution, diffusivity and energetics) as the longer
# format allows more pages.

using GLMakie, JLD2, GibbsSeaWater, ColorSchemes

cd(@__DIR__)
stable_output = "../dns_runs/stable_stepchange_nothing_600min/stable_analysis_and_field_snapshots.jld2"
lesscabbeling_output = "../dns_runs/lesscabbeling_stepchange_nothing_600min/lesscabbeling_analysis_and_field_snapshots.jld2"
cabbeling_output = "../dns_runs/cabbeling_stepchange_nothing_780min/cabbeling_analysis_and_field_snapshots.jld2"
isothermal_output = "../dns_runs/isothermal_stepchange_nothing_780min/isothermal_analysis_and_field_snapshots.jld2"
all_output = (isothermal_output, stable_output, lesscabbeling_output, cabbeling_output)
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
file = jldopen(cabbeling_output)

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

###########################################################################################
## New figures + replacing figure 3 with multiple figures
###########################################################################################

# Not exactly sure how where these will go so for now I will just make the figures

## S-Θ initial dns initial conditions.
# These should go in a table but maybe seeing them on S-Θ diagram will help. Much of this is
# the same as figure one so make sure that is run first.

haline_grad = get(ColorSchemes.haline, range(0, 1, length = 4))
fig = Figure(size = (1200, 500))
ax = Axis(fig[1, 1], xlabel = "S (gkg⁻¹)", ylabel = "Θ (°C)", title = "DNS initial conditions")

N = 1000
S_range, Θ_range = range(34.535, 34.705, length = N), range(-1.75, 0.75, length = N)
Θ_linear_initial = @. Θ_star + m_initial * (S_range - S_star)
S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
ρ = gsw_rho.(S_grid, Θ_grid, 0)
Θ_freezing = gsw_ct_freezing.(S_range, 0, 1)

S₀ᵘ = [34.69431424, 34.551, 34.568, 34.58]
Θ₀ᵘ = vcat(0.5, fill(-1.5, 3))
S_max = similar(S₀ᵘ)
Θ_max = similar(Θ₀ᵘ)
ρ_max_pred = Vector{Float64}(undef, 4)
for i ∈ eachindex(ρ_max)
    slope = (Θ₀ᵘ[i] - Θ_star) / (S₀ᵘ[i] - S_star)
    S_mix = range(S₀ᵘ[i], S_star, step = 0.000001)
    Θ_mix = @. Θ₀ᵘ[i] + (slope) * (S_mix - S₀ᵘ[i])
    ρ_mix = gsw_rho.(S_mix, Θ_mix, 0)
    ρ_max_pred[i], ρ_max_idx = findmax(ρ_mix)
    S_max[i], Θ_max[i] = S_mix[ρ_max_idx], Θ_mix[ρ_max_idx]
end
contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star],
         color = :grey, linewidth = 1,
         labelsize = 18, label = "Isopycnal at deep water")
lines!(ax, S_range, Θ_linear_initial,
        color = :grey, label = "Tangent at deep water", linestyle = :dashdot)
scatter!(ax, [S_star], [Θ_star], color = :red, label = "Deep water")
scatter!(ax, S₀ᵘ, Θ₀ᵘ, color = haline_grad, label = "Shallow water", marker = :rect)
# for i ∈ eachindex(S₀ᵘ)
#     lines!(ax, [S₀ᵘ[i], S_star], [Θ₀ᵘ[i], Θ_star], color = haline_grad[i], linewidth = 0.5,
#             linestyle = :dash, label = i == 3 ? "Mixing line" : nothing)
# end
# scatter!(ax, S_max[2:end], Θ_max[2:end], color = haline_grad[2:end], label = "Density maximum", marker = :utriangle)
# scatter!(ax, S_max[2:end], Θ_max[2:end], color = haline_grad[2:end], label = "Density maximum", marker = :cross)
# lines!(ax, S_range, Θ_freezing, label = "Freezing point", color = :blue)
Legend(fig[1, 2], ax)
fig

## Density evolution, perhaps pair this with figure 3 panel (a).
ρ_max = Vector{Float64}(undef, length(all_output))
ha_σ = Vector{Array}(undef, length(all_output))
zC = Vector{Array}(undef, length(all_output))
timestamps = Vector{Array}(undef, length(all_output))
for (i, file) ∈ enumerate(all_output)

    f = jldopen(file)
    timestamps[i] = f["dims/timestamps"]
    zC[i] = f["dims/z"]

    # Horizontally averaged profile
    ha_σ[i] = hcat([reshape(f["σ/σ_$(t)"], :) for t ∈ timestamps[i]]...)
    ρ_max[i] = f["attrib/ρ_max"]

    close(f)
end
ρ_max_model = [maximum(ha_σ[i][:, 2]) for i ∈ 1:4]
scatter_position = [1, 2, 3, 4]
expts = ["iso", "stable", "lesscab", "cab"]

fig = Figure(size = (500, 500))
ax = Axis(fig[1, 1], xlabel = "Experiment", xticks = (scatter_position, expts), ylabel = "σ₀ (kgm⁻³)")
scatter!(ax, scatter_position, ρ_max .- 1000, color = :steelblue, markersize = 15, label = "Theoretical predicted\nmaximum")
scatter!(ax, scatter_position, ρ_max_model .- 1000, color = :orange, marker = :xcross, markersize = 15, label = "HA maximum after\nmixing at intreface")
axislegend(ax, position = :lt)
fig

## Diffusivity, new panel for all depth integrated diffusivities but keep other colourmap
∫κₛ = Vector{Array}(undef, length(all_output))
for (i, file) ∈ enumerate(all_output)

    f = jldopen(file)
    ∫κₛ[i] = f["diffusivity"]["∫κₛ"]

    close(f)
end

restricted_time = 1:600
fig = Figure(size = (1000, 500))
expts = ["isothermal", "stable", "less cabbeling", "cabeling"]
ax = Axis(fig[1, 1], xlabel = "time (mins)", ylabel = "κ̅_eff (m²s⁻¹, log10)",
            title = "Depth integrated horizontally averaged salinity effective diffusivity",
            # xticks = time_ticks
            )
for (i, file) ∈ enumerate(∫κₛ)
    lines!(ax, timestamps[i][restricted_time] ./ 60, log10.(abs.(∫κₛ[i][restricted_time])),
            color = haline_grad[i], label = expts[i])
end
hlines!(ax, log10(1e-7), label = "Parameterised salinity diffusivity", linestyle = :dash,
        color = :black)
axislegend(ax, position = :rt, labelsize = 17)
xlims!(ax, (0, 600))
fig

## Energetics

file = jldopen(cabbeling_output)
ρ₀ = file["attrib/ρ₀"]
g = file["attrib/g"]
timestamps = file["dims/timestamps"]
 # the correction is due to incorrect scaling by area in the original computation
 bpe = file["energetics"]["∫Eb"]  * ρ₀ * (0.1^2 / 0.07^2)
  pe = file["energetics"]["∫Ep"]  * ρ₀
  ek = file["energetics"]["∫Ek"]  * ρ₀
   ε = file["energetics"]["∫ε"]   * ρ₀
∫gρw = file["energetics"]["∫gρw"] * ρ₀
  Δt = diff(timestamps)

ape = pe .- bpe
dₜek = diff(ek) ./ Δt
dₜpe = diff(pe) ./ Δt
dₜbpe = Φd = diff(bpe) ./ Δt
dₜape = diff(ape) ./ Δt
time_interp_mins = 0.5 * (timestamps[1:end-1] .+ timestamps[2:end]) / 60

close(file)

## (In) sanity check of energy budget
ε_interp = 0.5 * (ε[1:end-1] .+ ε[2:end])
∫gρw_interp = 0.5 * (∫gρw[1:end-1] .+ ∫gρw[2:end])
RHS = -ε_interp .- ∫gρw_interp
fig, ax = lines(time_interp, dₜek)
lines!(ax, time_interp, RHS)
fig

## I will still only show the energetics for the most active (i.e. cabbeling) case.
fig = Figure(size = (1000, 500))
ax = Axis(fig[1, 1], xlabel = "time (mins)", ylabel = "Joules", title = "(a) Energy reservoirs")
energy_scale = 1e5
lines!(ax, timestamps[restricted_time],  pe[restricted_time] .* energy_scale,
        label = "PE")
lines!(ax, timestamps[restricted_time], bpe[restricted_time] .* energy_scale,
        label = "BPE")
lines!(ax, timestamps[restricted_time], ape[restricted_time] .* energy_scale,
        label = "APE", color = :red)
axislegend(ax, position = :rb, orientation = :horizontal)
hidexdecorations!(ax, ticks = false, grid = false)
Label(fig[1, 1, Top()], halign = :left, L"\times 10^{-5}", fontsize = 16)
# time derivatives
ax2 = Axis(fig[2, 1], xlabel = "time (mins)", ylabel = "Watts",
          title = "(b) Time derivative of energy reservoirs")
dₜenergy_scale = 1e8
lines!(ax2, time_interp_mins[restricted_time],  dₜpe[restricted_time] .* dₜenergy_scale,
        label = "dₜPE")
lines!(ax2, time_interp_mins[restricted_time], dₜbpe[restricted_time] .* dₜenergy_scale,
        label = "dₜBPE")
lines!(ax2, time_interp_mins[restricted_time], dₜape[restricted_time] .* dₜenergy_scale,
        label = "dₜAPE", color = :red)
# lines!(ax3, time_interp_mins[1:400],  dₜek[1:400] .* energy_scale, label = "dₜEk", color = :green)
Label(fig[2, 1, Top()], halign = :left, L"\times 10^{-8}", fontsize = 16)
axislegend(ax2, position = :rt, orientation = :horizontal)

# topspinecolor = leftspinecolor = rightspinecolor = bottomspinecolor = :gray
# zoom_window = 1:100

# ax_inset = Axis(fig[2, 1],
#     width=Relative(0.68),
#     height=Relative(0.572),
#     halign=0.85,
#     valign=0.6,
#     backgroundcolor=:white,
#     xticklabelsize = 12,
#     yticklabelsize = 12;
#     topspinecolor, leftspinecolor, rightspinecolor, bottomspinecolor
#     )
# ylims!(ax_inset, -0.05, 0.4)
# lines!(ax_inset, time_interp_mins[zoom_window],  dₜpe[zoom_window] .* dₜenergy_scale)
# lines!(ax_inset, time_interp_mins[zoom_window], dₜbpe[zoom_window] .* dₜenergy_scale)
# lines!(ax_inset, time_interp_mins[zoom_window], dₜape[zoom_window] .* dₜenergy_scale, color = :red)
# hlines!(ax_inset, 0, linestyle = :dash, color = :grey)
# text!(ax2, 148, 3.85, text = L"\times 10^{-8}", fontsize = 12)
# translate!(ax_inset.scene, 0, 0, 10)
# # this needs separate translation as well, since it's drawn in the parent scene
# translate!(ax_inset.elements[:background], 0, 0, 9)

# using GLMakie.GeometryBasics
# ϵ = 1e-12
# rectangle_corners = Point2f[(0, minimum(dₜape)*energy_scale - ϵ),
#                             (time_interp_mins[zoom_window[end]], minimum(dₜape)*energy_scale-ϵ),
#                             (time_interp_mins[zoom_window[end]], maximum(dₜpe)*energy_scale+ϵ),
#                             (0, maximum(dₜpe)*energy_scale + ϵ)]
# poly!(ax2, rectangle_corners, strokecolor = :grey, strokewidth = 1, color = :transparent)
# xs = [time_interp_mins[310]]
# ys = [1]
# us = [time_interp_mins[zoom_window[end]+10] - time_interp_mins[315]]
# vs = [0]
# arrows!(ax2, xs, ys, us, vs, color = :gray)
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
