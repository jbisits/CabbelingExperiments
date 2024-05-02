using NCDatasets, JLD2, Statistics

tracers = "tracers.nc"

NCDataset(tracers) do ds

    time = ds[:time]

    for t ∈ eachindex(time)

        S = mean(ds[:S][:, :, :, t], dims = (1, 2))
        T = mean(ds[:T][:, :, :, t], dims = (1, 2))

        jldopen("cabbeling_profile.jld2") do file

            file["S_ha"] = reshape(S, :)
            file["T_ha"] = reshape(T, :)

        end

    end

end
