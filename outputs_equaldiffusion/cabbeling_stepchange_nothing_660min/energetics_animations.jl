using NCDatasets, CairoMakie, Printf, StatsBase, TwoLayerDirectNumericalShenanigans

"""
    function animate_reference_profile(computed_output::AbstractString; σ_binwidth = 0.0001)
Animate a pdf of the density as well as the reference profile for density
(density sorted from smallest to largest) `computed_output`.
"""
function animate_reference_profile(computed_output::AbstractString; σ_binwidth = 0.00015)

    NCDataset(computed_output) do ds

        t = ds["time"][:]

        pred_σₗ = ds.attrib["Predicted maximum density"]

        n = Observable(1)
        σ_extrema = [minimum(ds["σ"][:, :, :, 1]), pred_σₗ]
        σ_edges = σ_extrema[1]-σ_binwidth:σ_binwidth:σ_extrema[2]+σ_binwidth
        σ_hist = @lift fit(Histogram, reshape(ds["σ"][:, :, :, $n], :), σ_edges)
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (400, 800))
        ax = Axis(fig[1, 1], title = time_title)

        plot!(ax, σ_hist, color = :steelblue)
        ax.xlabel = "σ₀ (kgm⁻³)"
        ax.ylabel = "Frequency"
        vlines!(ax, pred_σₗ, color = :red, linestyle = :dash)

        ax2 = Axis(fig[2, 1], title = "Density reference profile",
                   xlabel = "σ₀ (kgm⁻³)", ylabel = "z* (m)")
        ax2.xticklabelrotation = π / 4
        σ_initial_reference_profile = sort(reshape(ds["σ"][:, :, :, 1], :), rev = true)
        z = range(-1, 0, length(σ_initial_reference_profile))
        lines!(ax2, σ_initial_reference_profile, z, color = :red, linestyle = :dot,
                label = "Initial reference profile")
        σ_reference_profile = @lift sort(reshape(ds["σ"][:, :, :, $n], :), rev = true)
        lines!(ax2, σ_reference_profile, z)
        vlines!(ax2, pred_σₗ, color = :red, linestyle = :dash, label = "Predicted σₗ")
        linkxaxes!(ax, ax2)
        hidexdecorations!(ax, grid = false)
        axislegend(ax2, position = :lb)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "reference_density.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing

end
"""
    function animate_reference_and_horizontal_average_profile(computed_output::AbstractString)
Animate the reference (sorted) density profile and the horizontally average profile.
"""
function animate_reference_and_horizontal_average_profile(computed_output::AbstractString)

    NCDataset(computed_output) do ds

        t = ds["time"][:]
        z_model = ds["zC"][:]

        pred_σₗ = ds.attrib["Predicted maximum density"]

        n = Observable(1)
        σₕ = @lift reshape(mean(ds["σ"][:, :, :, $n], dims = (1, 2)), :)
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (800, 400))
        ax = Axis(fig[1, 1], title = time_title)

        lines!(ax, σₕ, z_model, color = :steelblue)
        ax.subtitle = "Horizontally averaged profile"
        ax.xlabel = "σ₀ (kgm⁻³)"
        ax.ylabel = "z (m)"
        ax.xticklabelrotation = π / 4
        vlines!(ax, pred_σₗ, color = :red, linestyle = :dash)

        ax2 = Axis(fig[1, 2], title = "Density reference profile",
                   xlabel = "σ₀ (kgm⁻³)", ylabel = "z* (m)")
        ax2.xticklabelrotation = π / 4
        σ_initial_reference_profile = sort(reshape(ds["σ"][:, :, :, 1], :), rev = true)
        z = range(-1, 0, length(σ_initial_reference_profile))
        lines!(ax2, σ_initial_reference_profile, z, color = :red, linestyle = :dot,
                label = "Initial reference profile")
        σ_reference_profile = @lift sort(reshape(ds["σ"][:, :, :, $n], :), rev = true)
        lines!(ax2, σ_reference_profile, z)
        vlines!(ax2, pred_σₗ, color = :red, linestyle = :dash, label = "Predicted σₗ")
        linkyaxes!(ax, ax2)
        hideydecorations!(ax2, grid = false)
        axislegend(ax2, position = :lb)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "reference_and_horizontal_average_profile.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing

end
"""
   function ∫ρw!(energy_diagnostics::AbstractString, computed_output::AbstractString, velocities::AbstractString)
Volume integrated buoyancy flux.
"""
function ∫ρw!(energy_diagnostics::AbstractString, computed_output::AbstractString, velocities::AbstractString)

    NCDataset(computed_output) do ds

        var_key = "∫ρw"
        if !haskey(ds, var_key)

            time = ds[:time][:]
            g = -9.81
            dV = diff(ds[:xC][1:2]) .* diff(ds[:yC][1:2]) .* diff(ds[:zC][1:2])

            ∫ρw = similar(time)
            ds_velocities = NCDataset(velocities)
            for t ∈ eachindex(time)
                σ = ds[:σ][:, :, :, t]
                w = (ds_velocities[:w][:, :, 1:end-1, t] .+ ds_velocities[:w][:, :, 2:end, t]) / 2
                ∫ρw[t] = g * sum(σ .* w * dV[1])
            end

            NCDataset(energy_diagnostics, "a") do ds2

                defVar(ds2, var_key, ∫ρw, tuple("time"),
                    attrib = Dict("longname" => "Volume integrated buoyancy flux g∫ᵥρwdV."))

            end

            close(ds_velocities)
        end
    end

    return nothing
end
"""
    function dₜEb!(energy_diagnostics::AbstractString, computed_output::AbstractString)
Compute the time evolution of the background potential energy as in equation (12)
from Winters et al. (1995) and save to `energy_diagnostics`.
"""
function dₜEb!(energy_diagnostics::AbstractString, computed_output::AbstractString)

    NCDataset(computed_output) do ds

        var_key = "∫Eb"
        if !haskey(ds, var_key)

            time = ds[:time][:]
            g = -9.81
            x_length = length(ds[:xC])
            y_length = length(ds[:yC])
            z = repeat(ds[:zC][:], inner= x_length * y_length)
            dV = diff(ds[:xC][1:2]) .* diff(ds[:yC][1:2]) .* diff(ds[:zC][1:2])

            ∫Eb = similar(time)
            for t ∈ eachindex(time)
                σ = reshape(ds[:σ][:, :, :, t], :)
                p = sortperm(σ)
                ∫Eb[t] = g * sum(σ .* z[p] * dV[1])
            end

            NCDataset(energy_diagnostics, "a") do ds2

                defVar(ds2, var_key, ∫Eb, tuple("time"),
                    attrib = Dict("longname" => "Volume integrated background potential energy (∫ᵥgρz*dV)."))

            end

        end
    end

    return nothing

end
"""
    function dₜEp!(energy_diagnostics::AbstractString, computed_output::AbstractString)
Compute the time evolution of the potential energy  and save to `energy_diagnostics`.
"""
function dₜEp!(energy_diagnostics::AbstractString, computed_output::AbstractString)

    NCDataset(computed_output) do ds

        var_key = "∫Ep"
        if !haskey(ds, var_key)

            g = -9.81
            time = ds[:time][:]
            x_length = length(ds[:xC])
            y_length = length(ds[:yC])
            z_length = length(ds[:zC])
            z_grid = reshape(repeat(ds[:zC][:], inner= x_length * y_length),
                            (x_length, y_length, z_length))
            dV = diff(ds[:xC][1:2]) .* diff(ds[:yC][1:2]) .* diff(ds[:zC][1:2])
            ∫Ep = similar(time)

            for t ∈ eachindex(time)
                ∫Ep[t] = g * sum(ds[:σ][:, :, :, t] .* z_grid * dV[1])
            end

            NCDataset(energy_diagnostics, "a") do ds2

                defVar(ds2, var_key, ∫Ep, tuple("time"),
                    attrib = Dict("longname" => "Volume integrated potential energy (∫ᵥgρzdV)."))

            end

        end

    end

    return nothing

end
"""
    function Φᵢ!(energy_diagnostics::AbstractString, computed_output::AbstractString)
Compute the rate of conversion of internal energy to potential energy (Winters et al. (1995))
"""
    function Φᵢ!(energy_diagnostics::AbstractString, computed_output::AbstractString)

        NCDataset(computed_output) do ds

            var_key = "Φᵢ"
            if !haskey(ds, var_key)

                time = ds["time"][:]
                g = -9.81
                Δx = diff(ds["xC"][1:2])[1]
                Δy = diff(ds["yC"][1:2])[1]
                A = Δx * Δy * length(ds["xC"][:])^2
                σ = ds["σ"]
                κ = parse(Float64, ds.attrib["κₜ"][1:findfirst(' ', ds.attrib["κₜ"])])
                Φᵢ = similar(time)
                for t ∈ eachindex(time)
                    Φᵢ[t] = κ * g * A * (mean(σ[:, :, end, t]) -  mean(σ[:, :, 1, t]))
                end

                NCDataset(energy_diagnostics, "a") do ds2

                    defVar(ds2, var_key, Φᵢ, tuple("time"),
                        attrib = Dict("longname" => "Rate of conversion of internal energy to potential energy."))

                end

            end

        end

        return nothing

    end
"""
    function compute_energy_diagnostics!(energy_diagnostics::AbstractString,
                                         computed_output::AbstractString)
Save energy diagnostics to the file (must be a `.nc` file) `energy_diagnostics`.
"""
function compute_energy_diagnostics!(energy_diagnostics::AbstractString,
                                     computed_output::AbstractString,
                                     velocities::AbstractString)

    isfile(energy_diagnostics) ? nothing : makefile(energy_diagnostics, computed_output)

    dₜEb!(energy_diagnostics, computed_output)
    dₜEp!(energy_diagnostics, computed_output)
    ∫ρw!(energy_diagnostics, computed_output, velocities)
    Φᵢ!(energy_diagnostics, computed_output)

    return nothing

end
"""
    function makefile(filename::AbstractString, computed_output::AbstractString)
Create a file with `filename` that has the same dimensions as `computed_output`.
"""
function makefile(filename::AbstractString, computed_output::AbstractString)

    co = NCDataset(computed_output)

    NCDataset(joinpath(@__DIR__, filename), "c") do ds

        for key ∈ keys(co.dim)
            defDim(ds, key, length(co[key]))
            defVar(ds, co[key])
        end

        for key ∈ keys(co.attrib)
            ds.attrib[key] = co.attrib[key]
        end

        defVar(ds, co["∫ϵ"])
    end
    close(co)

    return nothing

end

velocities = "velocities.nc"
computed_output = "computed_output.nc"
# TLDNS.animate_density(computed_output, "σ")
# tracers = "tracers.nc"
# TLDNS.animate_tracers(tracers)
# animate_reference_profile(computed_output)
# energy_diagnostics = "energy_diagnostics.nc"
compute_energy_diagnostics!(energy_diagnostics, computed_output, velocities)
# animate_reference_and_horizontal_average_profile(computed_output)
