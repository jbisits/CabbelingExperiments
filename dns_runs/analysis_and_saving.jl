using NCDatasets, StatsBase, JLD2, GibbsSeaWater

"""
    function effective_diffusivity!(computed_output::AbstractString, tracers::AbstractString)
Compute the effective diffusivity for the salinity tracer `S` in `tracers`. The computed
values are then saved in `computed_output`.
"""
function effective_diffusivity!(computed_output::AbstractString, tracers::AbstractString)

    ds = NCDataset(tracers)

    time = ds[:time][:]
    Δt = diff(time)
    z = ds[:zC][:]
    Δz = diff(z)
    snapshots = eachindex(Δt)
    ∫κₛ = Vector{Float64}(undef, length(snapshots))
    κₛ_save = Array{Float64}(undef, length(z)-1, length(Δt))
    ∂S∂z_save = Array{Float64}(undef, length(z)-1, length(Δt))
    Fₛ = Array{Float64}(undef, length(z), length(Δt))

    for (i, t) ∈ enumerate(snapshots)

        S = [reshape(mean(ds[:S][:, :, :, t], dims = (1, 2)), :) reshape(mean(ds[:S][:, :, :, t+1], dims = (1, 2)), :)]
        sort!(S, dims = 1)
        ∫Sdz = cumsum(S * Δz, dims = 1)
        dₜ∫Sdz = diff(∫Sdz, dims = 2) ./ Δt[t]
        Fₛ[:, i] = dₜ∫Sdz
        ∂S∂z = diff(S, dims = 1) ./ Δz
        κₛ = dₜ∫Sdz[2:end] ./ ∂S∂z[:, 2]
        ∂S∂z_save[:, t] .= vec(∂S∂z[:, 2])
        replace!(κₛ, Inf => NaN)
        replace!(κₛ, 0 => NaN)
        replace!(κₛ, -Inf => NaN)
        find_nan = findall(.!isnan.(κₛ))
        κₛ_save[:, t] = κₛ
        ∫κₛ[i] = sum(κₛ[find_nan] .* Δz[find_nan]) / sum(Δz[find_nan])

    end

    close(ds)

    # Save to computed output
    NCDataset(computed_output, "a") do ds
        defVar()
    end

    return nothing
end
