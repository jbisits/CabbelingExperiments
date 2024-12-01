using NCDatasets, StatsBase, JLD2, GibbsSeaWater

tracers = "tracers.nc"
computed_output = "computed_output.nc"
velocities = "velocities.nc"

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

        defDim(ds, "diff_z", length(Δz))
        defDim(ds, "diff_t", length(Δt))

        defVar(ds, "Fₛ", κₛ, ("zC", diff_t),
                attrib = ("longname" => "Horizontally averaged vertical salt flux"))
        defVar(ds, "∂S∂z", ∂S∂z_save, (diff_z, diff_t),
                attrib = ("longname" => "Horizontally averaged vertical salt gradient"))
        defVar(ds, "κₛ", κₛ_save, (diff_z, diff_t),
                attrib = ("longname" => "Horizontally averaged effective diffusivity"))
        defVar(ds, "∫κₛ", ∫κₛ, (diff_t,),
                attrib = ("longname" => "Horizontally averaged depth integrated effective diffusivity"))
    end

    return nothing
end

"""
    function potential_and_background_potential_energy!(computed_output::AbstractString)
Compute and append the potential and background energy to `computed_output`. **Note** the
PE and BPE are both referenced to ``z = 0``, and the saved quantities are volume integrated.
"""
function potential_and_background_potential_energy!(computed_output::AbstractString)

    NCDataset(computed_output, "a") do ds

        t = ds[:time][:]
        SA = 0.1 * 0.1
        Δx = diff(ds[:xC][1:2])[1]
        Δy = diff(ds[:yC][1:2])[1]
        Δz = diff(ds[:zC][1:2])[1]
        ΔV = Δx * Δy * Δz

        V = cumsum(ones(length(reshape(ds[:σ][:, :, :, 1], :)))) * ΔV
        z✶ = V / SA

        find_num = findfirst('k', ds.attrib["Reference density"]) - 1
        ρ₀ = parse(Float64, ds.attrib["Reference density"][1:find_num])

        σ = ds[:σ]
        # background potential energy
        Eb = similar(t)
        g = 9.81
        for i ∈ eachindex(t)
            σᵢ = σ[:, :, :, i] .- ρ₀
            σᵢ_array = reshape(σᵢ, :)
            sort!(σᵢ_array, rev = true)
            Eb[i] = (g / ρ₀) * sum(σᵢ_array .* z✶ * ΔV)
        end

        # save
        defVar(ds, "∫Eb", Eb, (t,),
                attrib = ("longname" => "Volume integrated background potential energy"))

        # potential energy
        z = ds[:zC]
        z_ref0 = reverse(abs.(z))
        z_grid = reshape(repeat(z_ref0, inner = 124 * 124), (124, 124, 1400))
        Ep = similar(t)
        for i ∈ eachindex(t)
            σᵢ = σ[:, :, :, i] .- ρ₀
            Ep[i] = (g / ρ₀) * sum(σᵢ .* z_grid * ΔV)
        end

    end

    defVar(ds, "∫Ep", Ep, (t,),
            attrib = ("longname" => "Volume integrated background potential energy"))

    return nothing
end
"""
    function extract_and_save!(saved_data::AbstractString, computed_output::AbstractString)
Extract and save data needed for plots. The `saved_data` needs to be a `.jld2` file.
"""
function extract_and_save!(saved_data::AbstractString, computed_output::AbstractString, velocities::AbstractString)

    jldopen(saved_data, "w") do file

        NCDataset(computed_output) do ds

            file["∫Eb"] = ds["∫Eb"][:]
            file["∫Ep"] = ds["∫Ep"][:]
            file["∫Ek"] = ds["∫Eₖ"][:]
             file["∫ε"] = ds["∫ϵ"][:]
            file["∫κₛ"] = ds["∫κₛ"][:]
             file["κₛ"] = ds["κₛ"][:, :]

        end

    end

    return nothing
end


effective_diffusivity!(computed_output, tracers)
potential_and_background_potential_energy!(computed_output)

saved_data = "analysis_and_field_snapshots.jld2"
