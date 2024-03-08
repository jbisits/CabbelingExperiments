include("onedmodel.jl")
using .OneDModel

run_OneDModel(:cabbeling)

##
using TwoLayerDirectNumericalShenanigans, CairoMakie, JLD2, GibbsSeaWater

## Load output, cd = convective diffusivity value
output = joinpath(@__DIR__, "OneDModelOutput_cabbeling_cd1.jld2")
output = joinpath(@__DIR__, "OneDModelOutput_cabbeling_cd10.jld2")

S_ts, T_ts = FieldTimeSeries(output, "S"), FieldTimeSeries(output, "T")
σ₀_ts =  FieldTimeSeries(output, "σ")
time = S_ts.times
z = znodes(S_ts.grid, Center(), Center(), Center())
file = jldopen(output)
κ_ts = Array{Float64}(undef, length(z), length(time))
for (i, t) ∈ enumerate(string.(0:660))
    κ_ts[:, i] = file["timeseries"]["κ"][t][:]
end
close(file)

## Initial snapshots
visualise_snapshot(T_ts, "Θ (°C)", 1, 1, 1)
visualise_snapshot(S_ts, "S (gkg⁻ꜝ)", 1, 1, 1; colormap = :haline)

## Hovmoller plots
hplot = TShovmoller_plot(S_ts, T_ts)
save(joinpath(@__DIR__, "oned_TS_hovs_cd1.png"), hplot)

## Density
hplot_density = hovmoller_plot(σ₀_ts, "σ₀ (kgm⁻³)"; colormap = :dense)
save(joinpath(@__DIR__, "oned_density_hov_cd1.png"), hplot_density)
σ₀_max = maximum(σ₀_ts.data) # 1027.7103539732836, our predicted max density

## Diffusivity
fig, ax, hm = heatmap(time, z, κ_ts', colormap = :speed)
Colorbar(fig[1, 2], hm, label = "Diffusivity (m²s⁻¹)")
fig
save(joinpath(@__DIR__, "oned_diffusivity_hov_cd1.png"), fig)

## Fluxes
"""
    function tracer_flux(C_ts::FieldTimeSeries)
Calculate the `tracer_flux` from a `FieldTimeSeries`. The flux is calculated as the change
in tracer content within a fixed volume over time.
"""
function tracer_flux(C_ts::FieldTimeSeries)

    Δz = C_ts.grid.Δzᵃᵃᶜ
    time = C_ts.times
    Δt = diff(time)
    C_content = Array{Float64}(undef, length(z), length(time))

    for t ∈ eachindex(time)

        C_sorted = sort(reshape(C_ts[:, :, :, t], :))
        C_content[:, t] = cumsum(C_sorted * Δz)

    end

    flux = Array{Float64}(undef, length(z), length(time)-1)

    for t ∈ 1:length(time)-1
        flux[:, t] = diff(C_content[:, t:t+1], dims = 2) / Δt[1]
    end

    return reverse(flux, dims = 1)

end

Fₛ = tracer_flux(S_ts)
replace!(Fₛ, 0 => NaN)
fig, ax, hm = heatmap(time[2:end], z, Fₛ', colormap = :haline)
Colorbar(fig[1, 2], hm, label = "Salt flux ")
fig

Fₜ = tracer_flux(T_ts)
replace!(Fₜ, 0 => NaN)
fig, ax, hm = heatmap(time[2:end], z, Fₜ', colormap = :thermal)
Colorbar(fig[1, 2], hm, label = "Heat flux (Wm⁻²)")
fig

g = 9.81
α = gsw_alpha(34.6, -0.5, 0)
cₚ = 4000
β = gsw_beta(34.6, -0.5, 0)
J_b = @. -g * (α * Fₜ / cₚ - β * Fₛ)

fig, ax, hm = heatmap(time[2:end], z, J_b', colormap = :dense)
Colorbar(fig[1, 2], hm, label = "Buoyancy flux")
fig
