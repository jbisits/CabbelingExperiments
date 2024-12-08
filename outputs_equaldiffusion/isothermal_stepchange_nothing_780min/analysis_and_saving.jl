using NCDatasets, StatsBase, JLD2, CairoMakie
using TwoLayerDirectNumericalShenanigans: animate_density, animate_tracers

"""
    function effective_diffusivity!(computed_output::AbstractString, tracers::AbstractString)
Compute the effective diffusivity for the salinity tracer `S` in `tracers`. The computed
values are then saved in `computed_output`.
"""
function effective_diffusivity!(computed_output::AbstractString, tracers::AbstractString)

    ds = NCDataset(tracers)

    time = ds[:time][:]
    Δt = 0.5 * (time[1:end-1] .+ time[2:end])
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
        ∫Sdz = cumsum(S * Δz[1], dims = 1)
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

        haskey(ds.dim, "Δz") ? nothing : defDim(ds, "Δz", length(Δz))
        haskey(ds.dim, "Δt") ? nothing : defDim(ds, "Δt", length(Δt))

        Δt_vals = 0.5 * (ds[:time][1:end-1] .+ ds[:time][2:end])

          haskey(ds, "Δt") ? nothing : defVar(ds, "Δt", Δt_vals, ("Δt",))
          haskey(ds, "Fₛ") ? nothing : defVar(ds, "Fₛ", Fₛ, ("zC", "Δt"),
                                            attrib = ("longname" => "Horizontally averaged vertical salt flux"))
        haskey(ds, "∂S∂z") ? nothing : defVar(ds, "∂S∂z", ∂S∂z_save, ("Δz", "Δt"),
                                            attrib = ("longname" => "Horizontally averaged vertical salt gradient"))
          haskey(ds, "κₛ") ? nothing : defVar(ds, "κₛ", κₛ_save, ("Δz", "Δt"),
                                            attrib = ("longname" => "Horizontally averaged effective diffusivity"))
         haskey(ds, "∫κₛ") ? nothing : defVar(ds, "∫κₛ", ∫κₛ, ("Δt",),
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
        SA = 0.07 * 0.07
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

        # potential energy
        z = ds[:zC]
        Nx = ds.dim[:xC]
        Ny = ds.dim[:yC]
        Nz = ds.dim[:zC]
        z_ref0 = reverse(abs.(z))
        z_grid = reshape(repeat(z_ref0, inner = Nx * Ny), (Nx, Ny, Nz))
        Ep = similar(t)
        for i ∈ eachindex(t)
            σᵢ = σ[:, :, :, i] .- ρ₀
            Ep[i] = (g / ρ₀) * sum(σᵢ .* z_grid * ΔV)
        end


    # save
    haskey(ds, "∫Eb") ? nothing : defVar(ds, "∫Eb", Eb, ("time",),
                                        attrib = ("longname" => "Volume integrated background potential energy"))
    haskey(ds, "∫Ep") ? nothing : defVar(ds, "∫Ep", Ep, ("time",),
                                        attrib = ("longname" => "Volume integrated potential energy"))
    end

    return nothing
end
"""
    function buoyancy_flux!(computed_output::AbstractString, velocities::AbstractString)
Compute the buoyancy flux from model density and vertical velocity fields.
"""
function buoyancy_flux!(computed_output::AbstractString, velocities::AbstractString)

    ds_co = NCDataset(computed_output, "a")
    time = ds_co[:time][:]
    find_num = findfirst('k', ds_co.attrib["Reference density"]) - 1
    ρ₀ = parse(Float64, ds_co.attrib["Reference density"][1:find_num])
    ΔV = diff(ds_co[:xC][1:2])[1] * diff(ds_co[:yC][1:2])[1] * diff(ds_co[:zC][1:2])[1]
    ds_vel = NCDataset(velocities)

    g = 9.81
    ∫gρw = similar(time)
    for t ∈ eachindex(time)

        σ = ds_co[:σ][:, :, :, t]
        σ1 = @view σ[:, :, 1:end-1]
        σ2 = @view σ[:, :, 2:end]
        σ_interp = cat(σ[:, :, 1], 0.5 * (σ1 .+ σ2), σ[:, :, end], dims = 3)
        w = ds_vel[:w][:, :, :, t]

        ∫gρw[t] = (g / ρ₀) * sum(σ_interp .* w) * ΔV

    end
    haskey(ds_co, "∫gρw") ? nothing : defVar(ds_co, "∫gρw", ∫gρw, ("time",),
                                            attrib = Dict("longname" => "Volume integrated density flux in post processing"))
    close(ds_co)
    close(ds_vel)

    return nothing
end
"""
    function extract_and_save!(saved_data::AbstractString, computed_output::AbstractString,
                                velocities::AbstractString, snapshots::AbstractArray)
Extract and save data needed for plots. The `saved_data` needs to be a `.jld2` file.
**Note:** all the saved energetics data is per unit volume (i.e. I have divided by the
reference density).
"""
function extract_and_save!(saved_data::AbstractString, computed_output::AbstractString,
                            velocities::AbstractString, snapshots::AbstractArray,
                            interface_depths::AbstractArray)

    jldopen(saved_data, "w") do file

        NCDataset(computed_output) do ds

             file["dims/timestamps"] = ds["time"][:]
              file["dims/snapshots"] = ds["time"][snapshots]
                     file["dims/Δt"] = ds["Δt"][:]
                      file["dims/x"] = ds["xC"][:]
                      file["dims/y"] = ds["yC"][:]
                      file["dims/z"] = ds["zC"][:]

             file["energetics/∫Eb"] = ds["∫Eb"][:]
             file["energetics/∫Ep"] = ds["∫Ep"][:]
             file["energetics/∫Ek"] = ds["∫Eₖ"][:]
              file["energetics/∫ε"] = ds["∫ϵ"][:]
            file["energetics/∫gρw"] = ds["∫gρw"][:]

            file["diffusivity/∫κₛ"] = ds["∫κₛ"][:]
             file["diffusivity/κₛ"] = ds["κₛ"][:, :]

            find_num = findfirst('k', ds.attrib["Reference density"]) - 1
            ρ₀ = parse(Float64, ds.attrib["Reference density"][1:find_num])

            file["attrib/ρ₀"] = ρ₀
             file["attrib/g"] = 9.81

            t = ds["time"][:]
             for (i, t_) ∈ enumerate(t)
                file["σ/σ_$(t_)"] = mean(ds["σ"][:, :, :, i], dims = (1, 2))
            end

            for i ∈ snapshots
                file["σ/σ_xzslice/σ_$(t[i])"] = ds[:σ][:, 1, :, i]
                file["σ/σ_xyslice/σ_$(t[i])"] = ds[:σ][:, :, interface_depths[i], i]
            end
        end

        NCDataset(velocities) do ds

            file["dims/zF"] = ds["zF"][:]

            t = ds["time"][:]
            for i ∈ snapshots
                file["w/w_yzslice/w_$(t[i])"] = ds[:w][1, :, :, i] # 115 = end of y domain
                file["w/w_zmean/w_$(t[i])"] = mean(ds[:w][:, :, :, i], dims = 3)
            end

        end
    end

    return nothing
end

tracers = "tracers.nc"
computed_output = "computed_output.nc"
velocities = "velocities.nc"

@info "Computing diagnostics and saving to output"
effective_diffusivity!(computed_output, tracers)
potential_and_background_potential_energy!(computed_output)
buoyancy_flux!(computed_output, velocities)

saved_data = "isothermal_analysis_and_field_snapshots.jld2"
snapshots = 1:25
interface_depths = [821, 821, 821, 821, 821, 821, 821, 851, 841, 861,
                    861, 861, 871, 881, 891, 891, 901, 901, 901, 911,
                    920, 920, 925, 925, 925]

@info "Extracting and saving data required for plotting"
extract_and_save!(saved_data, computed_output, velocities, snapshots, interface_depths)

@info "Producing animations"
xslice = 57
yslice = 57
animate_density(computed_output, "σ"; xslice, yslice)
animate_tracers(tracers;  xslice, yslice)
