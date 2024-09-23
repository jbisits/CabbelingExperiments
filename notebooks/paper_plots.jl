# using CairoMakie, JLD2
using GLMakie, JLD2

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
z_xy_mid = z[Int(length(z)/2)] * ones(length(x), length(y))

snapshot = 10
slices = (density_xz_face1 = file["σ/σ_xzslice_face1/σ_$(60.0 * snapshot)"],
          density_xz_faceend = file["σ/σ_xzslice_faceend/σ_$(60.0 * snapshot)"],
          density_xz_mid = file["σ/σ_xzslice/σ_$(60.0 * snapshot)"],
          density_middle_slice = file["σ/σ_xyslice/σ_$(60.0 * snapshot)"],
          velocity_yz = file["w/w_yzslice/w_$(60.0 * snapshot)"][:, 1:end-1],
          velocity_yz_face1 = file["w/w_yzslice_face1/w_$(60.0 * snapshot)"][:, 1:end-1],
          velocity_yz_faceend = file["w/w_yzslice_faceend/w_$(60.0 * snapshot)"][:, 1:end-1],
          velocity_zmean = file["w/w_zmean/w_$(60.0 * snapshot)"][:, :, 1])

close(file)

## Make the figure
fig = Figure(size = (1600, 1600), fontsize = 18)

ax = Axis3(fig[1, 2],
           aspect=(1/3, 1/3, 1),
           xlabel = "x (m)",
           ylabel = "y (m)",
           zlabel = "z (m)",
           xlabeloffset = 50,
           ylabeloffset = 50,
           zlabeloffset = 75,
           limits = ((x[1], x[end]), (y[1], y[end]), (z[1], z[end])),
           elevation = π / 6.5,
           azimuth = 1.25π,
           xspinesvisible = false,
           zgridvisible = false,
           protrusions = 20,
           # persepectiveness = 0.1
           )

velocity_colorrange = maximum(abs.(extrema(slices.velocity_yz_face1))) * [-1, 1]
# velocity_colorrange = maximum(abs.(extrema(slices.velocity_zmean))) * [-1, 1]
sf_σ = surface!(ax, x_xz, y_xz_density, z_xz; color = slices.density_xz_face1, colormap = :dense)
sf_w = surface!(ax, x_xz_velocity, y_xz, z_yz; color = slices.velocity_yz_face1, colormap = :balance,
                colorrange = velocity_colorrange, ambient = (0.85, 0.85, 0.85), backlight = 1f0)
sf_σ = surface!(ax, x_xz, y_xz_density_end, z_xz; color = slices.density_xz_faceend, colormap = :dense)
sf_w = surface!(ax, x_xz_velocity_end, y_xz, z_yz; color = slices.velocity_yz_faceend, colormap = :balance,
                colorrange = velocity_colorrange, ambient = (0.85, 0.85, 0.85), backlight = 1f0)
surface!(ax, x, y, z_xy_top; color = slices.velocity_zmean, colormap = :balance, colorrange = velocity_colorrange)
surface!(ax, x, y, z_xy_mid; color = slices.density_middle_slice, colormap = :dense, overdraw = true) # add somthing to middle of domain

# Contours
# c_σ = contour!(ax, y, z, slices.density_xz; levels = 4, linewidth = 1, transformation = (:xz, y[1]),
#          transparency = true, color = :white, linestyle = :dash)


Colorbar(fig[1, 3], sf_σ, label = "σ₀ (kgm⁻³)")
Colorbar(fig[1, 1], sf_w, label = "w (ms⁻¹)", flipaxis = false)

ax.title = "DNS at t = $(times[snapshot + 1] / 60) minutes"
# fig[0, 1:2] = Label(fig, title; fontsize = 24, tellwidth = false, padding = (0, 0, -60, 0))

colgap!(fig.layout, 1, Relative(-0.1))
colgap!(fig.layout, 2, Relative(-0.1))

fig

##
save("dns_schematic.png", fig)

## multiple panels
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
z_xy_mid = z[Int(length(z)/2)] * ones(length(x), length(y))

snapshots = [0, 3, 6, 9]
slices = [(density_xz_face1 = file["σ/σ_xzslice_face1/σ_$(60.0 * snapshot)"],
          density_xz_faceend = file["σ/σ_xzslice_faceend/σ_$(60.0 * snapshot)"],
          density_xz_mid = file["σ/σ_xzslice/σ_$(60.0 * snapshot)"],
          density_middle_slice = file["σ/σ_xyslice/σ_$(60.0 * snapshot)"],
          velocity_yz = file["w/w_yzslice/w_$(60.0 * snapshot)"][:, 1:end-1],
          velocity_yz_face1 = file["w/w_yzslice_face1/w_$(60.0 * snapshot)"][:, 1:end-1],
          velocity_yz_faceend = file["w/w_yzslice_faceend/w_$(60.0 * snapshot)"][:, 1:end-1],
          velocity_zmean = file["w/w_zmean/w_$(60.0 * snapshot)"][:, :, 1]) for snapshot ∈ snapshots]

close(file)

##
fig = Figure(size = (1600, 1600), fontsize = 18)

ax = [Axis3(fig[i, j],
           aspect=(1/3, 1/3, 1),
           xlabel = "x (m)",
           ylabel = "y (m)",
           zlabel = "z (m)",
           xlabeloffset = 50,
           ylabeloffset = 50,
           zlabeloffset = 75,
           limits = ((x[1], x[end]), (y[1], y[end]), (z[1], z[end])),
           elevation = π / 6.5,
           azimuth = 1.25π,
           xspinesvisible = false,
           zgridvisible = false,
           protrusions = 20,
           # persepectiveness = 0.1
           ) for i ∈ 1:2, j ∈ 1:2]

velocity_colorrange = maximum([maximum(abs.(extrema(slices[i].velocity_yz_face1))) for i ∈ 1:4])  * [-1, 1]
density_colorrange = [minimum([minimum(slices[i].density_xz_face1) for i ∈ 1:4]),
                      maximum([maximum(slices[i].density_xz_face1) for i ∈ 1:4])]

for i ∈ eachindex(snapshots)
    sf_σ = surface!(ax[i], x_xz, y_xz_density, z_xz; color = slices[i].density_xz_face1,
                    colormap = :dense,
                    colorrange = density_colorrange
                    )
    sf_w = surface!(ax[i], x_xz_velocity, y_xz, z_yz; color = slices[i].velocity_yz_face1, colormap = :balance,
                    colorrange = velocity_colorrange,
                    ambient = (0.85, 0.85, 0.85), backlight = 1f0)
    sf_σ = surface!(ax[i], x_xz, y_xz_density_end, z_xz; color = slices[i].density_xz_faceend,
                    colormap = :dense,
                    colorrange = density_colorrange
                    )
    sf_w = surface!(ax[i], x_xz_velocity_end, y_xz, z_yz; color = slices[i].velocity_yz_faceend, colormap = :balance,
                    colorrange = velocity_colorrange,
                    ambient = (0.85, 0.85, 0.85), backlight = 1f0)
    surface!(ax[i], x, y, z_xy_top; color = slices[i].velocity_zmean, colormap = :balance,
            colorrange = velocity_colorrange
            )
    surface!(ax[i], x, y, z_xy_mid; color = slices[i].density_middle_slice, colormap = :dense,
           colorrange = density_colorrange,
            overdraw = true)
    ax[i].title = "t = $(times[snapshots[i] + 1] / 60) minutes"
    if i == 1
        Colorbar(fig[:, 3], sf_σ, label = "σ₀ (kgm⁻³)")
        Colorbar(fig[:, 0], sf_w, label = "w (ms⁻¹)", flipaxis = false)
    end
    if i > 2
        hidezdecorations!(ax[i], ticks = false)
    end
end
#colgap!(fig.layout, 1, Relative(-0.1))
colgap!(fig.layout, 2)
colgap!(fig.layout, 3, Relative(-0.1))
fig
##
save("dns_schematic_ts.png", fig)
