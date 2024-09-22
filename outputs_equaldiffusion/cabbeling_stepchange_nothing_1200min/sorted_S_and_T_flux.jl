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

    ∫Sdz = cumsum(S * Δz✶, dims = 1)
    dₜ∫Sdz = vec(diff(∫Sdz, dims = 2) / Δt[t])

    ∫Tdz = cumsum(T * Δz✶, dims = 1)
    dₜ∫Tdz = vec(diff(∫Tdz, dims = 2) / Δt[t])

    S1 = @view S[:, 1:end-1]
    S2 = @view S[:, 2:end]
    S_interp = 0.5 * (S1 .+ S2)
    T1 = @view T[:, 1:end-1]
    T2 = @view T[:, 2:end]
    T_interp = 0.5 * (T1 .+ T2)

    α = gsw_alpha.(S_interp, T_interp, 0)
    β = gsw_beta.(S_interp, T_interp, 0)

    βFₛ[t] = sum(β .* dₜ∫Sdz * Δz✶)
    αFₜ[t] = sum(α .* dₜ∫Tdz * Δz✶)

end
close(ds_tracers)
close(ds_computed_output)

jldopen(bflux, "a+") do file
    file["βFₛ_alt_2"] = βFₛ
    file["αFₜ_alt_2"] = αFₜ
end
