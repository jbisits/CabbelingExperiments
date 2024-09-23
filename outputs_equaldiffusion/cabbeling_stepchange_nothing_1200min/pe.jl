using NCDatasets, JLD2

computed_output = "computed_output.nc"
energetics_file = "cabbeling_energetics.jld2"
notebook_path = "/g/data/e14/jb2381/CabbelingExperiments/notebooks"
energetics_path = joinpath(notebook_path, energetics_file)

co_ds = NCDataset(computed_output)

Δx = diff(co_ds[:xC][1:2])[1]
Δy = diff(co_ds[:yC][1:2])[1]
Δz = diff(co_ds[:zC][1:2])[1]
ΔV = Δx * Δy * Δz

t = co_ds[:time][:]
z = co_ds[:zC][:]

find_num = findfirst('k', co_ds.attrib["Reference density"]) - 1
ρ₀ = parse(Float64, co_ds.attrib["Reference density"][1:find_num])
∫ϵ = co_ds[:∫ϵ][:]
∫Eₖ = co_ds[:∫Eₖ][:]

close(co_ds)

z_ref0 = reverse(abs.(z))
z_grid = reshape(repeat(z_ref0, inner = 124 * 124), (124, 124, 1400))
Ep = similar(t)
g = 9.81
for i ∈ eachindex(t)
    σᵢ = NCDataset(computed_output) do ds
            ds[:σ][:, :, :, i] .- ρ₀
    end
    Ep[i] = (g / ρ₀) * sum(σᵢ .* z_grid * ΔV)
end

jldopen(energetics_path, "a+") do file
    file["∫Ep_zref0_density_anomaly"] = Ep
    # file["ρ₀"] = ρ₀
    # file["∫ε"] = ∫ϵ
    # file["∫Ek"] = ∫Eₖ
end
