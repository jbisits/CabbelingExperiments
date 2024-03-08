include("onedmodel.jl")
using .OneDModel

run_OneDModel(:cabbeling)

##
using TwoLayerDirectNumericalShenanigans, CairoMakie, JLD2

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

    for t ∈ 1:leanght(time)-1
        flux[:, t] = diff(C_content[:, t:t+1]) / Δt
    end

    return flux

end

Fₛ = tracer_flux(S_ts)
