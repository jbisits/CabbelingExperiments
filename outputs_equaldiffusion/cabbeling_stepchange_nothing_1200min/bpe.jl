using NCDatasets, JLD2

computed_output = "computed_output.nc"
energetics_file = "cabbeling_energetics.jld2"
notebook_path = "/g/data/e14/jb2381/CabbelingExperiments/notebooks"
energetics_path = joinpath(notebook_path, energetics_file)

co_ds = NCDataset(computed_output)

t = co_ds[:time][:]
SA = 0.1 * 0.1
Δx = diff(co_ds[:xC][1:2])[1]
Δy = diff(co_ds[:yC][1:2])[1]
Δz = diff(co_ds[:zC][1:2])[1]
ΔV = Δx * Δy * Δz
ΔA = Δx * Δy

V = cumsum(ones(length(reshape(co_ds[:σ][:, :, :, 1], :)))) * ΔV
z✶ = V / SA

find_num = findfirst('k', co_ds.attrib["Reference density"]) - 1
ρ₀ = parse(Float64, co_ds.attrib["Reference density"][1:find_num])

close(co_ds)

Eb = similar(t)
g = 9.81
for i ∈ eachindex(t)
    σᵢ = NCDataset(computed_output) do ds
            ds[:σ][:, :, :, i] .- ρ₀
    end
    σᵢ_array = reshape(σᵢ, :)
    sort!(σᵢ_array, rev = true)
    Eb[i] = (g / ρ₀) * sum(σᵢ_array .* z✶ * ΔV)
end
jldopen(energetics_path, "a+") do file
    file["∫Eb′z✶"] = Eb
end
