using CairoMakie, NCDatasets, Printf

isothermal = "outputs_equaldiffusion/isothermal_stepchange_nothing_660min/computed_output.nc"
cabbeling = "outputs_equaldiffusion/cabbeling_stepchange_nothing_660min/computed_output.nc"

function animate_density(isothermal::AbstractString, cabbeling::AbstractString, variable::AbstractString="σ";
                         yslice = 52)

    NCDataset(isothermal) do iso

        x_iso = iso["xC"][:]
        z_iso = iso["zC"][:]
        NCDataset(cabbeling) do cab

            x_cab = cab["xC"][:]
            z_cab = cab["zC"][:]
            t = cab["time"][:]

            pred_max_density = cab.attrib["Predicted maximum density"]

            n = Observable(1)
            σ_iso = @lift iso[variable][:, yslice, :, $n]
            σ_cab = @lift cab[variable][:, yslice, :, $n]
            time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

            fig = Figure(size = (1000, 500))
            ax = [Axis(fig[1, i], title = i == 1 ? time_title : "") for i ∈ 1:2]

            colormap = cgrad(:dense)[2:end-1]
            colorrange = (minimum(cab[variable][:, :, :, 1]), pred_max_density)
            lowclip = cgrad(:dense)[1]
            highclip = cgrad(:dense)[end]

            # isothermal
            heatmap!(ax[1], x_iso, z_iso, σ_iso; colorrange, colormap, lowclip, highclip)
            ax[1].xlabel = "x (m)"
            ax[1].ylabel = "z (m)"

            # cabbeling
            hm = heatmap!(ax[2], x_cab, z_cab, σ_cab; colorrange, colormap, lowclip, highclip)
            ax[2].xlabel = "x (m)"
            ax[2].ylabel = "z (m)"
            Colorbar(fig[1, 3], hm, label = "σ₀ (kgm⁻³)")

            linkyaxes!(ax[1], ax[2])
            hideydecorations!(ax[2], ticks = false)

            frames = eachindex(t)
            record(fig, joinpath(pwd(), "cab_iso_density.mp4"),
            frames, framerate=8) do i
                msg = string("Plotting frame ", i, " of ", frames[end])
                print(msg * " \r")
                n[] = i
            end

        end

    end

    return nothing
end

animate_density(isothermal, cabbeling)
