using NCDatasets, JLD2, GibbsSeaWater

tracers = "tracers.nc"
computed_output = "computed_output.nc"
bflux = "buoyancy_flux_interp_face.jld2"

ds_tracers = NCDataset(tracers)
ds_computed_output = NCDataset(computed_output)

time = ds_tracers[:time][:]
Δt = diff(t)
SA = 0.1 * 0.1
Δx = diff(ds_computed_output[:xC][1:2])[1]
Δy = diff(ds_computed_output[:yC][1:2])[1]
Δz = diff(ds_computed_output[:zC][1:2])[1]
ΔV = Δx * Δy * Δz
ΔA = Δx * Δy

V = cumsum(ones(length(reshape(ds_computed_output[:σ][:, :, :, 1], :)))) * ΔV
z✶ = V / SA
Δz✶ = diff(z✶)

find_num = findfirst('k', ds_computed_output.attrib["Reference density"]) - 1
ρ₀ = parse(Float64, ds_computed_output.attrib["Reference density"][1:find_num])
g = 9.81

Fₛ = similar(Δt)
βFₛ = similar(Δt)
Fₜ = similar(Δt)
αFₜ = similar(Δt)

for t ∈ eachindex(Δt)

    S = reshape(ds_tracers[:S][:, :, :, t:t+1], :)
    sort!(S, dims = 1, rev = true)
    T = reshape(ds_tracers[:T][:, :, :, t:t+1], :)
    sort!(T, dims = 1, rev = true)
    α = gsw_alpha.(S, T, 0)
    β = gsw_beta.(S, T, 0)

    ## Salt without β
    ∫Sdz = cumsum(S * Δz✶, dims = 1)
    dₜ∫Sdz = vec(diff(∫Sdz, dims = 2) / Δt[t])
    Fₛ[t] = g * sum(dₜ∫Sdz .* z✶ * ΔV)

    ## Overwrite for memory efficieny
    ∫Sdz = cumsum(β .* S * Δz✶, dims = 1)
    dₜ∫Sdz = vec(diff(∫Sdz, dims = 2) / Δt[t])
    βFₛ[t] = g * sum(∫Sdz .* z✶ * ΔV)

    ## Temperature withour α
    ∫Tdz = cumsum(T * Δz✶, dims = 1)
    dₜ∫Tdz = vec(diff(∫Tdz, dims = 2) / Δ[t])
    Fₜ[t] = g * sum(dₜ∫Tdz .* z✶ * ΔV)

    ## Overwrite for memory efficiency
    ∫Tdz = cumsum(α .* T * Δz✶, dims = 1)
    dₜ∫Tdz = vec(diff(∫Tdz, dims = 2) / Δ[t])
    αFₜ[t] = g * sum(dₜ∫Tdz .* z✶ * ΔV)

end
close(ds_tracers)
close(ds_computed_output)

jldopen(bflux, "a+") do file
    file["Fₛ"] = Fₛ
    file["Fₜ"] = Fₜ
    file["βFₛ"] = βFₛ
    file["αFₜ"] = αFₜ
    file["ρ₀"] = ρ₀
end
