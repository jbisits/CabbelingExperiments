using NCDatasets, CairoMakie, Printf, StatsBase

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

"""
function ρw!(energy_diagnostics::AbstractString, computed_output::AbstractString)

    NCDataset(computed_output) do ds

        time = ds[:time][1:2]
        dz = diff(ds[:zC][1:2])
        dA = diff(ds[:xC][1:2]) .* diff(ds[:yC][1:2])

        ρw = similar(time)
        for t ∈ eachindex(time)
            σ = sum(ds[:σ][:, :, 350:1050, t] * dz[1], dims = 3) # sum over density range [-0.75, -0.25].
            ρw[t] = sum(σ * dA[1])
        end
        ρw = ρw .- ρw[1]

        NCDataset(energy_diagnostics, "a") do ds2

            defVar(ds2, "ρw", ρw, tuple("time"),
                   attrib = Dict("longname" => "Vertical buoyancy flux integrated from z = [-1, -0.25]."))

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

        time = ds[:time][1:2]
        g = -9.81
        x_length = length(ds[:xC])
        y_length = length(ds[:yC])
        z = repeat(ds[:zC][:], inner= x_length * y_length)
        dV = diff(ds[:xC][1:2]) .* diff(ds[:yC][1:2]) .* diff(ds[:zC][1:2])

        ∫Eb = similar(time)
        for t ∈ eachindex(time)
            σ = reshape(ds[:σ][:, :, :, t], :)
            p = sortperm(σ)
            ∫Eb[t] = g * sum(σ[p] .* z[p] * dV[1])
        end

        NCDataset(energy_diagnostics, "a") do ds2

            defVar(ds2, "∫Eb", ∫Eb, tuple("time"),
                   attrib = Dict("longname" => "Volume integrated background potential energy (∫ᵥgρz*dV)."))

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

        g = -9.81
        time = ds[:time][1:2]
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

            defVar(ds2, "∫Ep", ∫Ep, tuple("time"),
                   attrib = Dict("longname" => "Volume integrated potential energy (∫ᵥgρzdV)."))

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
                                     computed_output::AbstractString)

    isfile(energy_diagnostics) ? nothing : makefile(energy_diagnostics, computed_output)

    dₜEb!(energy_diagnostics, computed_output)
    dₜEp!(energy_diagnostics, computed_output)
    ρw!(energy_diagnostics, computed_output)

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

computed_output = "computed_output.nc"
energy_diagnostics = "energy_diagnostics.nc"
animate_reference_profile(computed_output)
compute_energy_diagnostics!(energy_diagnostics, computed_output)
