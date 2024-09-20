using NCDatasets, JLD2, GibbsSeaWater

tracers = "tracers.nc"
computed_output = "computed_output.nc"
bflux = "buoyancy_flux_interp_face.jld2"

ds_tracers = NCDataset(tracers)
ds_computed_output = NCDataset(computed_output)

time = ds_tracers[:time][:]
Δt = diff(time)
SA = 0.1 * 0.1
Δx = diff(ds_computed_output[:xC][1:2])[1]
Δy = diff(ds_computed_output[:yC][1:2])[1]
Δz = diff(ds_computed_output[:zC][1:2])[1]
ΔV = Δx * Δy * Δz
ΔA = Δx * Δy

V = cumsum(ones(length(reshape(ds_computed_output[:σ][:, :, :, 1], :)))) * ΔV
z✶ = V / SA
Δz✶ = diff(z✶)[1]

find_num = findfirst('k', ds_computed_output.attrib["Reference density"]) - 1
ρ₀ = parse(Float64, ds_computed_output.attrib["Reference density"][1:find_num])
g = 9.81

Fₛ = similar(Δt)
βFₛ = similar(Δt)
Fₜ = similar(Δt)
αFₜ = similar(Δt)

for t ∈ eachindex(Δt)

    S = [reshape(ds_tracers[:S][:, :, :, t], :) reshape(ds_tracers[:S][:, :, :, t+1], :)]
    sort!(S, dims = 1, rev = true)
    T = [reshape(ds_tracers[:T][:, :, :, t], :) reshape(ds_tracers[:T][:, :, :, t+1], :)]
    sort!(T, dims = 1)
    α = gsw_alpha.(S, T, 0)
    β = gsw_beta.(S, T, 0)

    ∫Sdz = cumsum(β .* S * Δz✶, dims = 1)
    dₜ∫Sdz = vec(diff(∫Sdz, dims = 2) / Δt[t])
    βFₛ[t] = sum(dₜ∫Sdz * Δz✶)

    ∫Tdz = cumsum(α .* T * Δz✶, dims = 1)
    dₜ∫Tdz = vec(diff(∫Tdz, dims = 2) / Δt[t])
    αFₜ[t] = sum(dₜ∫Tdz * Δz✶)

end
close(ds_tracers)
close(ds_computed_output)

jldopen(bflux, "a+") do file
    file["βFₛ"] = βFₛ
    file["αFₜ"] = αFₜ
end
