using GLMakie, JLD2, GibbsSeaWater, ColorSchemes, StatsBase
using SeawaterPolynomials: TEOS10EquationOfState, total_density, haline_contraction, thermal_expansion
using TwoLayerDirectNumericalShenanigans: REFERENCE_DENSITY
eos = TEOS10EquationOfState(reference_density = REFERENCE_DENSITY)

# Load output
cd(@__DIR__)
stable_output = "../dns_runs/stable_stepchange_nothing_600min/stable_analysis_and_field_snapshots.jld2"
lesscabbeling_output = "../dns_runs/lesscabbeling_stepchange_nothing_600min/lesscabbeling_analysis_and_field_snapshots.jld2"
cabbeling_output = "../dns_runs/cabbeling_stepchange_nothing_780min/cabbeling_analysis_and_field_snapshots.jld2"
isothermal_output = "../dns_runs/isothermal_stepchange_nothing_780min/isothermal_analysis_and_field_snapshots.jld2"
unstable_output = "../dns_runs/unstable_stepchange_nothing_600min/unstable_analysis_and_field_snapshots.jld2"
all_output = (isothermal_output, stable_output, lesscabbeling_output, cabbeling_output, unstable_output)

## Figure theme
haline_grad = get(ColorSchemes.haline, range(0, 5 / 5.3, length = length(all_output)))
rmn_label_stability = ["I, stable to cabbeling", "II, stable to cabbeling", "III, unstable to cabbeling",
                        "IV, unstable to cabbeling", "V, statically unstable"]
rmn_label = ["I", "II", "III", "IV", "V"]
markers = [:utriangle, :dtriangle, :cross, :cross, :rect]
markersize = 19
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
new_theme = merge(theme_latexfonts(), publication_theme)
set_theme!(new_theme)

## Figure 1: schematic
S✶, Θ✶ = 34.7, 0.5
S₀ᵘ = 34.58
Θ₀ᵘ = -1.5
slope = (Θ₀ᵘ - Θ✶) / (S₀ᵘ - S✶)
S_mix = range(S₀ᵘ, S✶, step = 0.000001)
Θ_mix = @. Θ₀ᵘ + (slope) * (S_mix - S₀ᵘ)
ρ_mix = gsw_rho.(S_mix, Θ_mix, 0)
ρ_max, ρ_max_idx = findmax(ρ_mix)
S_max, Θ_max = S_mix[ρ_max_idx], Θ_mix[ρ_max_idx]
Δρ_mix = ρ_max - gsw_rho(S✶, Θ✶, 0)

N = 1000
plotting_offset = 10
S_range, Θ_range = range(34.55, 34.705, length = N), range(-2, 1, length = N)
S_grid, Θ_grid = ones(N-(plotting_offset - 1)) .* S_range[plotting_offset:end]',
                 ones(N-(plotting_offset - 1))' .* Θ_range[plotting_offset:end]
ρ = gsw_rho.(S_grid, Θ_grid, 0)
ρ✶ = gsw_rho(S✶, Θ✶, 0)
ρ_s = gsw_rho(S₀ᵘ, Θ₀ᵘ, 0)
find_Θ = findfirst(Θ_range .> -1.5)
find_S = findfirst(ρ[find_Θ, :] .> ρ✶)
S_iso, Θ_iso = S_range[find_S], Θ_range[find_Θ]
gsw_rho(S_iso, Θ_iso, 0)
αₗ, βₗ = gsw_alpha(S✶, Θ✶, 0), gsw_beta(S✶, Θ✶, 0)
m_initial = βₗ / αₗ
Θ_linear_initial = @. Θ✶ + m_initial * (S_range - S✶)
αₘ, βₘ = gsw_alpha(S_max, Θ_max, 0), gsw_beta(S_max, Θ_max, 0)
m = βₘ / αₘ
Θ_linear = @. Θ_max + m * (S_range - S_max)

# Vertical profile
Nz = 200
z = range(-10, 0, length = 2*Nz)

Sˡ, Θˡ = fill(S✶, Nz), fill(Θ✶, Nz)
σ₀ˡ = gsw_sigma0(Sˡ[1], Θˡ[1])
Sᵘ, Θᵘ = fill(S₀ᵘ, Nz), fill(Θ₀ᵘ, Nz)
σ₀ᵘ = gsw_sigma0(Sᵘ[1], Θᵘ[1])
S, Θ = vcat(Sˡ, Sᵘ), vcat(Θˡ, Θᵘ)
σ₀ = gsw_sigma0.(S, Θ)

Θ_mix_profile = (Θᵘ[1] + Θˡ[1]) / 2
S_mix_profile = (Sᵘ[1] + Sˡ[1]) / 2
mixing_interface = 190:210
S[mixing_interface] = fill(S_mix_profile, length(mixing_interface))
Θ[mixing_interface] = fill(Θ_mix_profile, length(mixing_interface))
σ₀_mix = gsw_sigma0.(S, Θ)
Δσ₀ᵘ = σ₀_mix[mixing_interface[10]] - σ₀ᵘ
Δσ₀ˡ = σ₀_mix[mixing_interface[10]] - σ₀ˡ
zplot = range(-10, 10, length = 2*Nz)
upper_σ_mix = @. σ₀ᵘ + (Δσ₀ᵘ / 2) *( 1 + tanh(400* ((z - (-5)) / 10)))
lower_σ_mix = @. σ₀ˡ + (Δσ₀ˡ / 2) *( 1 + tanh(100* ((z - (-5)) / 10)))
σ_mix = vcat(lower_σ_mix[10:Nz+5], reverse(upper_σ_mix[1:Nz+4]))

fig = Figure(size = (1500, 600), px_per_unit=16)

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
lines!(ax, σ_mix[Nz-23:Nz+4], z[Nz-23:Nz+4], label = "Mixed line", color = :purple, linewidth = 2)
lines!(ax, fill(σ₀ˡ, Nz-20), z[1:Nz-20], label = "Deep water", color = :red, linewidth = 2)
lines!(ax, fill(σ₀ᵘ, Nz-2), z[Nz+3:end], label = "Shallow water", color = :blue, linewidth = 2)
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
text!(ax, maximum(σ₀_mix)-0.0008, z[mixing_interface[end]+10], align = (:right, :center),
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
contour!(ax2, S_range[plotting_offset:end], Θ_range[plotting_offset:end], ρ'; levels = [ρ_s, ρ✶, ρ_max],
         color = [:blue, :red, :magenta], linestyle = :dot, linewidth = 2,
         labelsize = 18, label = "Isopycnals")
lines!(ax2, S_mix, Θ_mix, color = :purple, label = "Mixed water", linewidth = 0.8)
scatter!(ax2, [S✶], [Θ✶]; color = :red, label = "Deep water", markersize)
scatter!(ax2, [S₀ᵘ], [Θ₀ᵘ]; color = :blue, label = "Shallow water", markersize)
lines!(ax2, S_range[plotting_offset:end], Θ_linear_initial[plotting_offset:end],
        color = (:red, 0.5), label = "Tangent at\ndeep water", linestyle = :dashdot)
scatter!(ax2, S_max, Θ_max; color = :magenta, label = "Maximum\ndensity", markersize)

# wedge fill in
ρ✶_isopycnal = findall(ρ' .≈ ρ✶)
S_idx = [ρ✶_isopycnal[i][1] for i ∈ eachindex(ρ✶_isopycnal)]
Θ_idx = [ρ✶_isopycnal[i][2] for i ∈ eachindex(ρ✶_isopycnal)]
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
text!(ax2, S✶-0.002, Θ✶, align = (:right, :bottom), text = "Deep water", color = :red)
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

## Figure 2: DNS flow evolution
# Load data to create figure
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
           titlefont = :regular,
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
Label(fig[0, :], "Flow evolution in Experiment IV", fontsize = 22, font = :bold)
rowgap!(fig.layout, 1, Relative(0.04))
rowgap!(fig.layout, 2, Relative(0.08))
fig
##
save("dns_schematic_ts_horizontal.png", fig)

## Figure 3 S-Θ initial conditions
fig = Figure(size = (1000, 500))
ax = Axis(fig[1, 1], xlabel = "S (gkg⁻¹)", ylabel = "Θ (°C)", title = "Experiment initial conditions")

N = 1000
S_range, Θ_range = range(34.535, 34.705, length = N), range(-1.75, 0.75, length = N)
Θ_linear_initial = @. Θ✶ + m_initial * (S_range - S✶)
S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
ρ = gsw_rho.(S_grid, Θ_grid, 0)
ρ✶ = gsw_rho(S✶, Θ✶, 0)
find_Θ = findfirst(Θ_range .> -1.5)
find_S = findfirst(ρ[find_Θ, :] .> ρ✶) - 1
S_iso = S_range[find_S]
Θ_freezing = gsw_ct_freezing.(S_range, 0, 1)

S₀ᵘ = [34.69431424, 34.551, 34.568, 34.58, 34.59]
Θ₀ᵘ = vcat(0.5, fill(-1.5, length(S₀ᵘ)-1))
S_max = similar(S₀ᵘ)
Θ_max = similar(Θ₀ᵘ)
ρ_max_pred = Vector{Float64}(undef, length(Θ_max))

for i ∈ eachindex(ρ_max)
    slope = (Θ₀ᵘ[i] - Θ✶) / (S₀ᵘ[i] - S✶)
    S_mix = range(S₀ᵘ[i], S✶, step = 0.000001)
    Θ_mix = @. Θ₀ᵘ[i] + (slope) * (S_mix - S₀ᵘ[i])
    ρ_mix = gsw_rho.(S_mix, Θ_mix, 0)
    ρ_max_pred[i], ρ_max_idx = findmax(ρ_mix)
    S_max[i], Θ_max[i] = S_mix[ρ_max_idx], Θ_mix[ρ_max_idx]
end

contour!(ax, S_range, Θ_range, ρ'; levels = [ρ✶],
         color = :grey, linewidth = 1,
         labelsize = 18, label = "Isopycnal at deep water")
lines!(ax, S_range, Θ_linear_initial,
        color = :grey, label = "Tangent at deep water", linestyle = :dashdot)
scatter!(ax, [S✶], [Θ✶]; color = :red, label = "Deep water", markersize)
for i ∈ eachindex(all_output)
    scatter!(ax, S₀ᵘ[i], Θ₀ᵘ[i]; color = i, colormap = :haline, colorrange = (0, 5.3),
            label = rmn_label_stability[i], marker = markers[i], markersize)
end
axislegend(ax, position = :lt)
fig
##
save("ics.png", fig)

## Figure 4: maximum density + density evolution
timestamps = Vector{Array}(undef, length(all_output))
zC = Vector{Array}(undef, length(all_output))
ha_σ = Vector{Array}(undef, length(all_output))
ρ_max = Vector{Float64}(undef, length(all_output))
for (i, file) ∈ enumerate(all_output)

    f = jldopen(file)
    timestamps[i] = f["dims/timestamps"]
    zC[i] = f["dims/z"]

    # Horizontally averaged profile
    ha_σ[i] = hcat([reshape(f["σ/σ_$(t)"], :) for t ∈ timestamps[i]]...)
    ρ_max[i] = f["attrib/ρ_max"]

    close(f)
end
ρ₀ = jldopen(all_output[1]) do f
        f["attrib/ρ₀"]
end
N = 100
ρ_max_theoretical = Vector{Float64}(undef, N)
Δρ_plot = similar(ρ_max_theoretical)
model_ρ✶ = total_density(Θ✶, S✶, 0, eos)
for (i, _S₀ᵘ) ∈ enumerate(range(34.52, 34.6, length = N))
    slope = (Θ₀ᵘ[end] - Θ✶) / (_S₀ᵘ - S✶)
    S_mix = range(_S₀ᵘ, S✶, step = 0.000001)
    Θ_mix = @. Θ₀ᵘ[end] + (slope) * (S_mix - _S₀ᵘ)
    eos_vec = fill(eos, length(S_mix))
    zero_vec = zeros(length(S_mix))
    ρ_mix = total_density.(Θ_mix, S_mix, zero_vec, eos_vec)
    ρ_max_theoretical[i] = maximum(ρ_mix)
    Δρ_plot[i] = total_density(Θ₀ᵘ[end], _S₀ᵘ, 0, eos) - model_ρ✶
end

ρ_max_model = [maximum(ha_σ[i][:, 2]) for i ∈ eachindex(all_output)]

Δρ = Vector{Float64}(undef, length(all_output))
for i ∈ eachindex(Δρ)
    Δρ[i] = total_density(Θ₀ᵘ[i], S₀ᵘ[i], 0, eos) - model_ρ✶
end
ΔS_xaxis = range(34.52, 34.6, length = N) .- S✶
ΔS_ics = S₀ᵘ .- S✶
snapshots = [1, 50, 100, 400, 600] .+ 1
expt = 4
σ_anomaly = 1000
fig = Figure(size = (1100, 500))
ax = Axis(fig[1, 1],
            title = "(a) Maximum density",
            xlabel = L"Initial $\Delta S$ (gkg$^{-1}$)",
            ylabel = L"$σ_{0}'$ (kgm$^{-3}$)",
            )
hlines!(ax, model_ρ✶ .- σ_anomaly, color = :grey, label = "Initial deep water density")
lines!(ax, ΔS_xaxis[30:end], ρ_max_theoretical[30:end] .- σ_anomaly,
        label = L"$ρ_{\mathrm{max}}$ prediction", color = :red, alpha = 0.5)
for i ∈ 2:length(Δρ)
    scatter!(ax, ΔS_ics[i], ρ_max_model[i] .- σ_anomaly; color = i, colormap = :haline,
            colorrange = (1, 5), label = rmn_label_stability[i],
            marker = markers[i], markersize)
end
axislegend(ax, position = :lt)

ax2 = Axis(fig[1, 2], xlabel = L"$σ_{0}'$ (kgm$^{-3}$)", ylabel = "z (m)")
vlines!(ax2, ρ_max[expt] .-  σ_anomaly,
color = :red, linestyle = :dash, label = L"$\rho_{\mathrm{max}}$ prediction")
lines!(ax2, ha_σ[expt][:, 1] .-  σ_anomaly, zC[expt],
label = "Initial profile", color = :black)
for (i, t) ∈ enumerate(snapshots)
    ax2.title = "(b) Experiment $(rmn_label[expt])"
    lines!(ax2, ha_σ[expt][:, t] .- σ_anomaly, zC[expt],
            label = "$(timestamps[expt][t] / 60) mins", alpha = 0.75)
end
axislegend(ax2, position = :lb)
fig
##
save("density_evolution.png", fig)

## Figure 5: effective diffusivity
∫κₛ = Vector{Array}(undef, length(all_output))
κₛ = Array{Float64}(undef, 1650, 781)
zC = Vector{Array}(undef, length(all_output))
timestamps = Vector{Array}(undef, length(all_output))
for (i, file) ∈ enumerate(all_output)

    f = jldopen(file)

    timestamps[i] = f["dims/timestamps"]
    zC[i] = f["dims/z"]
    ∫κₛ[i] = f["diffusivity"]["∫κₛ"]
    i == 4 ? κₛ = reverse(f["diffusivity/κₛ"], dims = 1) : nothing

    close(f)
end

restricted_time = 1:600
fig = Figure(size = (1000, 1000))
ax = Axis(fig[1, 1], xlabel = "time (mins)", ylabel = "z(m)",
        title = "(a) Experiment IV horizontally averaged salinity effective diffusivity")
hm = heatmap!(ax, timestamps[1][restricted_time] ./ 60, zC[4], log10.(abs.(κₛ[:, restricted_time]')),
                colorrange = (log10(1e-8), log10(1)), colormap = :tempo)
Colorbar(fig[1, 2], hm, label = L"$κ_{eff}$ (m²s$^{-1}$, log10)")
hidexdecorations!(ax, ticks = false, grid = false)
expts = ["isothermal", "stable", "less cabbeling", "cabbeling", "unstable"]
ax2 = Axis(fig[2, 1], xlabel = "time (mins)", ylabel = L"$\overline{κ_{eff}}$ (m²s$^{-1}$, log10)",
            title = "(b) Depth integrated horizontally averaged salinity effective diffusivity",
            )
for (i, file) ∈ enumerate(∫κₛ)
    lines!(ax2, timestamps[i][restricted_time] ./ 60, log10.(abs.(∫κₛ[i][restricted_time])),
            color = haline_grad[i], label = rmn_label_stability[i])
end
hlines!(ax2, log10(1e-7), label = "Parameterised salinity diffusivity", linestyle = :dash,
        color = :black)
axislegend(ax2, position = :rt, labelsize = 17)
xlims!(ax2, (0, 600))
linkxaxes!(ax, ax2)
fig
##
save("diffusivity.png", fig)

## Energetics
pe = Vector{Array}(undef, length(all_output))
bpe = similar(pe)
ape = similar(pe)
dₜpe = similar(pe)
dₜbpe = similar(pe)
dₜape = similar(pe)
dₜek = similar(pe)
ε = similar(pe)
∫gρw = similar(pe)
timestamps = similar(pe)
time_interp_mins = similar(pe)
for (i, file) ∈ enumerate(all_output)

    jldopen(file) do f

    ρ₀ = f["attrib/ρ₀"]
    g = f["attrib/g"]
    T = f["dims/timestamps"][end]
    timestamps[i] = f["dims/timestamps"]
    # the correction is due to incorrect scaling by area in the original computation
    pe₀ = f["energetics"]["∫Ep"][1]  * ρ₀
    scale_correction = i == 4 ? (0.1^2 / 0.07^2) : 1
    if i < 5
         pe[i] = (pe₀ .- f["energetics"]["∫Ep"] * ρ₀) ./ pe₀
        bpe[i] = (pe₀ .- f["energetics"]["∫Eb"] * scale_correction * ρ₀) ./ pe₀
    else
         pe[i] = (f["energetics"]["∫Ep"] * ρ₀ .- pe₀) ./ pe₀
        bpe[i] = (f["energetics"]["∫Eb"] * ρ₀ .- pe₀) ./ pe₀
    end
    ek = f["energetics"]["∫Ek"]
    ε[i] = f["energetics"]["∫ε"]
    ∫gρw[i] = f["energetics"]["∫gρw"]
    Δt = diff(timestamps[i])

    ape[i] = pe[i] .- bpe[i]
    dₜek[i] = diff(ek) ./ Δt
    dₜpe[i] = diff(pe[i]) ./ Δt
    dₜbpe[i] = Φd = diff(bpe[i]) ./ Δt
    dₜape[i] = diff(ape[i]) ./ Δt
    time_interp_mins[i] = 0.5 * (timestamps[i][1:end-1] .+ timestamps[i][2:end]) / 60

    end
end

## Figure A1: Batchelor sacle for expt 5
# Batchelor scale for unstable
Ba_unstable = Array{Float64}(undef, length(timestamps[5]))
jldopen(unstable_output) do file
    Ba_unstable .= file["attrib/Ba_array"] * 1e3
end
find_under_Ba = findall(Ba_unstable .< 0.6)
fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1], title = "Batchelor length in Experiment V",
            xlabel = "time (mins)", ylabel = "Ba (mm)")
hlines!(ax, 0.6, color = :red, linestyle = :dash, label = "Model resolution")
lines!(ax, timestamps[5][1:50] ./ 60, Ba_unstable[1:50], label = "Minimum local Batchelor scale in space")
scatter!(ax, timestamps[5][find_under_Ba] ./ 60, Ba_unstable[find_under_Ba], color = :orange,
         markersize = 8, label = "Batachelor resolution not resolved")
ylims!(ax, 0, 1)
axislegend(ax, position = :rb)
fig
##
save("batch_length_V.png", fig)

## Figure A2: energy budget
restricted_time = 1:600
fig = Figure(size = (900, 1800))
ax = [Axis(fig[i, 1], ylabel = "") for i ∈ 1:5]
for i ∈ 1:5
    ε_interp = 0.5 * (ε[i][1:end-1] .+ ε[i][2:end])
    ∫gρw_interp = 0.5 * (∫gρw[i][1:end-1] .+ ∫gρw[i][2:end])
    LHS = dₜek[i][restricted_time]
    RHS = -ε_interp .- ∫gρw_interp
    lines!(ax[i], LHS, label = L"\frac{1}{\rho_{0}}\frac{\mathrm{d}}{\mathrm{d}t}KE", linewidth = 2)
    lines!(ax[i], RHS[restricted_time], label = L"∫bw\mathrm{d}V - \frac{ε}{\rho_{0}}", linestyle = :dash, linewidth = 2)
    MAE = mean(abs.(LHS .- RHS[restricted_time]))
    ax[i].title = "Expt $(rmn_label[i]) - MAE = $(round(MAE, digits = 15))"
    if i ∈ 1:2
        ylims!(ax[i], -1e-14, 1e-14)
    end
    if i < 5
        hidexdecorations!(ax[i], grid = false, ticks = false)
    else
        ax[i].xlabel = "time (min)"
    end
end
Legend(fig[6, :], ax[1], orientation = :horizontal)
fig
##
save("energy_budgets.png", fig)

## Figure 6: cabbeling APE
expts = (4)
restricted_time = 1:600
fig = Figure(size = (1000, 500))
ax = Axis(fig[1, 1], xlabel = "time (mins)", ylabel = "Non-dimensional energy",
        title = "(a) Energy reservoirs")
linestyles = [:solid, :dash]
for (i, expt) ∈ enumerate(expts)
    if i < 2
    lines!(ax, timestamps[expt][restricted_time] / 60,  pe[expt][restricted_time],
            label = L"\mathcal{PE}", color = Makie.wong_colors()[1], linestyle = linestyles[i])
    lines!(ax, timestamps[expt][restricted_time] / 60, bpe[expt][restricted_time],
            label = L"\mathcal{BPE}", color = Makie.wong_colors()[2], linestyle = linestyles[i])
    end
    lines!(ax, timestamps[expt][restricted_time] / 60, ape[expt][restricted_time];
            label = L"\mathcal{APE}", color = :red, linestyle = linestyles[i])
end
axislegend(ax, position = :rc, orientation = :horizontal)
hidexdecorations!(ax, ticks = false, grid = false)

# time derivatives
ax2 = Axis(fig[2, 1], xlabel = "time (mins)", ylabel = L"Rate $(s^{-1})$",
          title = "(b) Time derivative of energy reservoirs")
for (i, expt) ∈ enumerate(expts)
    if i < 2
    lines!(ax2, time_interp_mins[expt][restricted_time],  dₜpe[expt][restricted_time],
            label = L"\mathrm{d}_{t}\mathcal{PE}", color = Makie.wong_colors()[1], linestyle = linestyles[i])
    lines!(ax2, time_interp_mins[expt][restricted_time], dₜbpe[expt][restricted_time],
            label = L"\mathrm{d}_{t}\mathcal{BPE}", color = Makie.wong_colors()[2], linestyle = linestyles[i])
    end
    lines!(ax2, time_interp_mins[expt][restricted_time], dₜape[expt][restricted_time],
            label = L"\mathrm{d}_{t}\mathcal{APE}", color = :red, linestyle = linestyles[i])
end
axislegend(ax2, position = :rt, orientation = :horizontal)

linkxaxes!(ax, ax2)

fig
##
save("cabbeling_ape.png", fig)

## Figure 7: close up of energy and diffusivity for IV
expt = 4
restricted_time = 250:400
∫gρw_interp = 0.5 * (∫gρw[expt][1:end-1] .+ ∫gρw[expt][2:end])
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "time (mins)", ylabel = L"Rate $(s^{-1})$",
          title = "Experiment IV, t = $(restricted_time[1])-$(restricted_time[end])mins",
          yaxisposition = :left,
          rightspinevisible = false)
l1 = lines!(ax, time_interp_mins[expt][restricted_time], dₜape[expt][restricted_time],
            color = :red)
l2 = lines!(ax, time_interp_mins[expt][restricted_time], dₜpe[expt][restricted_time],
            color = Makie.wong_colors()[1])
l3 = lines!(ax, time_interp_mins[expt][restricted_time], dₜbpe[expt][restricted_time],
            color = Makie.wong_colors()[2])
ax2 = Axis(fig[1, 1], ylabel = L"$\overline{κ_{eff}}$ (m²s$^{-1}$, log10)",
            yaxisposition = :right,
            yticklabelcolor = haline_grad[expt],
            rightspinecolor = haline_grad[expt],
            ytickcolor = haline_grad[expt],
            ylabelcolor = haline_grad[expt],
            leftspinevisible = false)
l4 = lines!(ax2, time_interp_mins[expt][restricted_time], log10.(abs.(∫κₛ[expt][restricted_time])),
            color = haline_grad[expt], label = rmn_label_stability[expt])
Legend(fig[2, 1], [l1, l2, l3, l4],
        [L"\mathrm{d}_{t}\mathcal{APE}", L"\mathrm{d}_{t}\mathcal{PE}",
        L"\mathrm{d}_{t}\mathcal{BPE}", L"$\overline{κ_{eff}}$"], orientation = :horizontal)
fig
##
save("ape_effdiff_IV.png", fig)

## Figure 9: density difference vs mixing rate
fig = Figure(size = (900, 700))
ylimits = (log10(1e-8), log10(1e-2))
ax = Axis(fig[1, 1], title = "(a) Static density difference vs mixing rate",
            xlabel = L"Intial $\Delta \rho$ (kgm$^{-3}$)",
            ylabel = L"$\langle\overline{κ_{eff}}\rangle_{t}$ (m²s$^{-1}$, log10)",
            limits = ((-0.035, 0.0089), ylimits),
            )
band_unstable = range(1e-5, 0.0099, length = 20)
vlines!(ax, 0, label = "Static stability threshold", color = :blue)
band!(band_unstable, ylimits..., color = (:blue, 0.25), label = "Statically unstable")
for i ∈ 1:length(Δρ)
    scatter!(ax, Δρ[i], log10(mean_∫κₛ[i]); color = haline_grad[i],
            marker = markers[i], markersize)
end
axislegend(ax, position = :lt)

ax2 = Axis(fig[2, 1], title = "(b) Cabbeling density difference vs mixing rate",
            xlabel = L"Intial $\Delta\rho'$ (kgm$^{-3}$)",
            ylabel = L"$\langle\overline{κ_{eff}}\rangle_{t}$ (m²s$^{-1}$, log10)",
            limits = ((nothing, 0.011) , ylimits)
            )
vlines!(ax2, 0, color= :red, label = "Stable to cabbeling")
xband = range(0, 0.011, length = 20)
bd = band!(ax2, xband, ylimits..., color = (:red, 0.25), label = "Unstable to cabbeling")
axislegend(ax2, position = :rb)#, backgroundcolor = :white)
sc = []
for i ∈ 1:length(Δρ)
    scatter!(ax2, Δρ_mix_upper[i], log10(mean_∫κₛ[i]); color = haline_grad[i],
            marker = markers[i], markersize, label = rmn_label_stability[i])
end
marker_label = [MarkerElement(color = haline_grad[i], marker = markers[i]; markersize) for i ∈ 1:5]
Legend(fig[3, :], marker_label, rmn_label_stability, orientation = :horizontal, nbanks = 2)
fig
##
save("eff_diff_instability_2panel.png", fig)

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
##

## Non-dimensional parameters
S₀ᵘ = [34.69431424, 34.551, 34.568, 34.58, 34.59]
Θ₀ᵘ = vcat(0.5, fill(-1.5, length(S₀ᵘ)-1))
function compute_R_ρ(S₀ᵘ, Θ₀ᵘ, eos; S₀ˡ = 34.7, Θ₀ˡ = 0.5)

    S_u = S_g = S₀ᵘ
    S_l = S_f = S₀ˡ
    T_u = T_f = Θ₀ᵘ
    T_l = T_g = Θ₀ˡ

    ρ_u = total_density(T_u, S_u, -0.5, eos)
    ρ_l = total_density(T_l, S_l, -0.5, eos)
    ρ_f = total_density(T_f, S_f, -0.5, eos)
    ρ_g = total_density(T_g, S_g, -0.5, eos)

    R_ρ = (0.5 * (ρ_f - ρ_u) + 0.5 * (ρ_l - ρ_g)) / (0.5 * (ρ_f - ρ_l) + 0.5 * (ρ_u - ρ_g))

    return R_ρ
end
for i ∈ eachindex(all_output)
    println(round(compute_R_ρ(S₀ᵘ[i], Θ₀ᵘ[i], eos), digits = 2))
end
function compute_Ra(S₀ᵘ, Θ₀ᵘ, eos;
                      S₀ˡ = 34.7, Θ₀ˡ = 0.5, g = 9.81, L = 0.5, κ = 1e-7, ν = 1e-6)

    S̄ = 0.5 * (S₀ˡ - S₀ᵘ)
    T̄ = 0.5 * (Θ₀ˡ - Θ₀ᵘ)
    ΔS = S₀ˡ - S₀ᵘ
    ΔT = Θ₀ˡ - Θ₀ᵘ
    S_u = S_g = S₀ᵘ
    S_l = S_f = S₀ˡ
    T_u = T_f = Θ₀ᵘ
    T_l = T_g = Θ₀ˡ

    ρ_u = total_density(T_u, S_u, -0.5, eos)
    ρ_l = total_density(T_l, S_l, -0.5, eos)
    ρ_f = total_density(T_f, S_f, -0.5, eos)
    ρ_g = total_density(T_g, S_g, -0.5, eos)

    ρ_m = total_density(T̄, S̄, -0.5, eos)

    # McDougall (1981)
    α̃ = (0.5 * (ρ_f - ρ_g) - 0.5 * (ρ_l - ρ_u)) / (ρ_m * ΔT)
    β̃ = (0.5 * (ρ_f - ρ_g) + 0.5 * (ρ_l - ρ_u)) / (ρ_m * ΔS)

    # α̃ = thermal_expansion(T̄, S̄, -0.5, eos)
    # β̃ = haline_contraction(T̄, S̄, -0.5, eos)

    Ra_S = (g * β̃ * ΔS * L^3) / (ν * κ)
    Ra_T = (g * α̃ * ΔT * L^3) / (ν * κ)

    return Ra_S, Ra_T
end
for i ∈ eachindex(all_output)
    println(round.(compute_Ra(S₀ᵘ[i], Θ₀ᵘ[i], eos)))
end
rayleigh_numbers = [compute_Ra(S₀ᵘ[i], Θ₀ᵘ[i], eos) for i ∈ eachindex(S₀ᵘ)]
