using NCDatasets, JLD2, GibbsSeaWater, Statistics

"""
    function mean_and_reshape_C(C1, C2, α1, α2)
Take the horizontal average of Cᵢ .* αᵢ, srt then return matrix. Here `C` is a tracer and
`α` is an expansion coefficient (haline or thermal depending on tracer).
"""
function mean_and_reshape_C(C1, C2, α1, α2)

    Cα1 = mean(C1 .* α1, dims = (1, 2))
    Cα2 = mean(C2 .* α2, dims = (1, 2))

    return [reshape(Cα1, :) reshape(Cα2, :)]
end
"""
    function compute_flux(C, Δz, Δt)
Compute the tracer flux between the two columns of the `C` matrix. This is done by
dₜ∫Cdz where ∫Cdz is the cumulative integral (of tracer content) to each depth.
"""
function compute_flux(C, Δz, Δt)

    sort!(C, dims = 1)
    ∫Cdz = cumsum(C * Δz, dims = 1)
    dₜ∫Cdz = diff(∫Cdz, dims = 2) ./ Δt

    return dₜ∫Cdz

end

cab_flux_file = "cabbeling_fluxes_and_diff_longer_run.jld2"
notebook_path = "/g/data/e14/jb2381/CabbelingExperiments/notebooks"
cab_flux_path = joinpath(notebook_path, cab_flux_file)

tracers = "tracers.nc"
ds = NCDataset(tracers)

time = ds[:time][:]
Δt = diff(time)
z = ds[:zC][:]
Δz = diff(z[1:2])[1]
snapshots = eachindex(Δt)
Fₛ = Array{Float64}(undef, length(z), length(Δt))
Fₜ = Array{Float64}(undef, length(z), length(Δt))

for (i, t) ∈ enumerate(snapshots)

    S1 = ds[:S][:, :, :, t]
    S2 = ds[:S][:, :, :, t+1]
    T1 = ds[:T][:, :, :, t]
    T2 = ds[:T][:, :, :, t+1]

    ## Salinity
    α1 = gsw_beta.(S1, T1, 0)
    α2 = gsw_beta.(S2, T2, 0)
    C = mean_and_reshape_C(S1, S2, α1, α2)
    Fₛ[:, i] = compute_flux(C, Δz, Δt[t])

    ## Temperature
    α1 = gsw_alpha.(S1, T1, 0)
    α2 = gsw_alpha.(S2, T2, 0)
    C = mean_and_reshape_C(T1, T2, α1, α2)
    Fₜ[:, i] = compute_flux(C, Δz, Δt[t])

end

close(ds)

jldopen(cab_flux_path, "a+") do file
    file["Fₛ_with_β"] = Fₛ
    file["Fₜ_with_α"] = Fₜ
end
