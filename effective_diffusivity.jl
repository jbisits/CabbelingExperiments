using NCDatasets, JLD2

tracers = "outputs_equaldiffusion/" .* ["isothermal_stepchangelineargradient_nothing_300min",
                                        "isothermal_stepchange_nothing_660min",
                                        "cabbeling_stepchangelineargradient_nothing_300min",
                                        "cabbeling_stepchange_nothing_660min"] .* "/tracers.nc"

save_names = ["isothermal/salinity_gradient", "isothermal/no_gradient",
              "cabbeling/salinity_gradient", "cabbeling/no_gradient"]

effective_diffusivity = "effective_diffusivity.jld2"
jldopen(effective_diffusivity, "w") do file
    file["time/salinity_gradient"] = NCDataset(tracers[1]) do ds
        ds[:time][:]
    end
    file["time/no_gradient"] = NCDataset(tracers[2]) do ds
        ds[:time][:]
    end
end

for (i, tracer) ∈ enumerate(tracers)

    ds = NCDataset(tracer)

    time = ds[:time][:]
    Δt = diff(time)
    ΔV = diff(ds[:xC][1:2])[1] * diff(ds[:yC][1:2])[1] * diff(ds[:zC][1:2])[1]
    V = (1:length(reshape(ds[:S][:, :, :, 1], :))) * ΔV
    SA = 0.1 * 0.1
    ẑ = V / SA
    Δẑ = diff(ẑ[1:2])[1]
    ∫κₛ = similar(Δt)

    for t ∈ eachindex(Δt)
        S = [reshape(ds[:S][:, :, :, t], :) reshape(ds[:S][:, :, :, t+1], :)]
        sort!(S, dims = 1)
        ∫Sdẑ = cumsum(S * Δẑ, dims = 1)
        dₜ∫Sdẑ = diff(∫Sdẑ, dims = 2) ./ Δt[t]
        ∂S∂ẑ = diff(S, dims = 1) ./ Δẑ
        κₛ = dₜ∫Sdẑ[2:end] ./ ∂S∂ẑ[:, 2]
        replace!(κₛ, Inf => 0)
        replace!(κₛ, NaN => 0)
        replace!(κₛ, -Inf => 0)
        ∫κₛ[t] = sum(κₛ * Δẑ)
    end

    close(ds)

    jldopen(effective_diffusivity, "a+") do file
        file[save_names[i]] = ∫κₛ[t]
    end

end
