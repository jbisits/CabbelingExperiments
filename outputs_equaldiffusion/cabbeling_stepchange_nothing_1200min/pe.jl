using NCDatasets, JLD2

computed_output = "computed_output.nc"
cab_flux_file = "cabbeling_fluxes_and_diff_longer_run.jld2"
notebook_path = "/g/data/e14/jb2381/CabbelingExperiments/notebooks"
cab_flux_path = joinpath(notebook_path, cab_flux_file)

co_ds = NCDataset(computed_output)

Δx = diff(co_ds[:xC][1:2])[1]
Δy = diff(co_ds[:yC][1:2])[1]
Δz = diff(co_ds[:zC][1:2])[1]
ΔV = Δx * Δy * Δz

t = co_ds[:time][:]
z = co_ds[:zC][:]

find_num = findfirst('k', co_ds.attrib["Reference density"]) - 1
ρ₀ = parse(Float64, co_ds.attrib["Reference density"][1:find_num])

close(co_ds)

z_grid = reshape(repeat(z, inner= 124 * 124), (124, 124, 1400))
Ep = similar(t)
g = 9.81
for i ∈ eachindex(t)
    σᵢ = NCDataset(computed_output) do ds
            ds[:σ][:, :, :, i]
    end
    Ep[i] = (g / ρ₀) * sum(σᵢ .* z_grid * ΔV)
end
jldopen(cab_flux_path, "a+") do file
    file["∫Ep_alternate"] = Ep
end
