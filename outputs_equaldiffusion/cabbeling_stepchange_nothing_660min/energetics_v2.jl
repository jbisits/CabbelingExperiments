using NCDatasets, JLD2

cab_tracers = "tracers.nc"
cab_computed_output = "computed_output.nc"
cab_energetics = "energy_diagnostics.nc"

co_ds = NCDataset(cab_computed_output)
t = co_ds[:time][:]
find_num = findfirst('k', co_ds.attrib["Reference density"]) - 1
ρ₀ = parse(Float64, co_ds.attrib["Reference density"][1:find_num])
dV = diff(co_ds[:xC][1:2]) .* diff(co_ds[:yC][1:2]) .* diff(co_ds[:zC][1:2])
z = co_ds[:zC][:]
close(co_ds)

en_ds = NCDataset(cab_energetics)
∫ϵ = ρ₀ .* en_ds[:∫ϵ][:]
∫Eₖ = ρ₀ .* en_ds[:∫Eₖ][:]
close(en_ds)

z_grid = reshape(repeat(z, inner= 124 * 124), (124, 124, 1400))
z_ = repeat(z, inner = 124*124)

Ep = Vector{Float64}(undef, length(t))
Eb = similar(Ep)
for i ∈ eachindex(t)
    σᵢ = NCDataset(cab_computed_output) do ds
            ds[:σ][:, :, :, i]
    end
    σᵢ_array = reshape(σᵢ, :)
    Ep[i] = 9.81 * sum(σᵢ .* z_grid * dV[1])
    p = sortperm(σᵢ_array, rev = true) # Missing the `rev = true` in the full calculation!!!!
    Eb[i] = 9.81 * sum(σᵢ_array .* z_[p] * dV[1])
end

cabbeling_energy = "cabbeling_energetics.jld2"
jldopen(cabbeling_energy, "w") do file
    file["Eb"] = Eb
    file["Ep"] = Ep
    file["Ek"] = ∫Eₖ
    file["ϵ"] = ∫ϵ
    file["time"] = t
end
