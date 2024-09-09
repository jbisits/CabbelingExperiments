using NCDatasets, JLD2

computed_output = "computed_output.nc"
cab_flux_file = "cabbeling_fluxes_and_diff_longer_run.jld2"
notebook_path = "/g/data/e14/jb2381/CabbelingExperiments/notebooks"
cab_flux_path = joinpath(notebook_path, cab_flux_file)

co_ds = NCDataset(computed_output)
ΔV = diff(co_ds[:xC][1:2])[1] * diff(co_ds[:yC][1:2])[1] * diff(co_ds[:zC][1:2])[1]
t = co_ds[:time][:]
z = co_ds[:zC][:]
close(co_ds)
z_ = repeat(z, inner = 124*124)

Eb = similar(t)
for i ∈ eachindex(t)
    σᵢ = NCDataset(computed_output) do ds
            ds[:σ][:, :, :, i]
    end
    σᵢ_array = reshape(σᵢ, :)
    p = sortperm(σᵢ_array, rev = true) # Missing the `rev = true` in the full calculation!!!!
    Eb[i] = 9.81 * sum(σᵢ_array .* z_[p] * ΔV)
end

jldopen(cab_flux_path, "a+") do file
    file["∫Eb"] = ∫Eb
end
