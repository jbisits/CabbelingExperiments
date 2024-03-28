using NCDatasets

"""
    function salt_checks!(file::AbstractString, tracers::AbstractString)
Compute and save some checks for the background potential energy of the salinity field and
the salt flux.
"""
function salt_checks!(file::AbstractString, tracers::AbstractString)

    ds = NCDataset(tracers)
    time = ds[:time][:]
    x_length = length(ds[:xC])
    y_length = length(ds[:yC])
    z = repeat(ds[:zC][:], inner= x_length * y_length)

    Δt = diff(time)
    Δx = diff(ds[:xC][1:2])[1]
    Δy = diff(ds[:yC][1:2])[1]
    Δz = diff(ds[:zC][1:2])[1]

    ΔV = Δx * Δy * abs(Δz)

    C = ds[:S]
    C_size = size(C)
    C_length = prod(C_size[1:3])
    V = (1:C_length) * ΔV
    SA = 0.1 * 0.1 # Surface area, would be better if not hard coded
    equivalent_z = -V / SA

    NCDataset(file, "c") do new_ds

        defDim(new_ds, "cumulative_volume", length(V))
        defDim(new_ds, "time", length(time))
        defDim(new_ds, "time_derivative", length(Δt))
        defVar(new_ds, "time", time, tuple("time"))
        defVar(new_ds, "time_derivative", cumsum(Δt), tuple("time_derivative"))
        defVar(new_ds, "volume", V, tuple("cumulative_volume"))
        defVar(new_ds, "equivalent_z", equivalent_z, tuple("cumulative_volume"))
        defVar(new_ds, "∫Sz✶dV", Float64, tuple("time"),
                attrib = Dict("longname" => "Background salinity potential energy (∫Sz✶dV)."))
        defVar(new_ds, "∫SdV", Float64, tuple("volume", "time"),
                attrib = Dict("longname" => "Cumulative salt content ∫SdV"))
        defVar(new_ds, "dₜ∫SdV", Float64, tuple("volume", "time_derivative"),
                attrib = Dict("longname" => "Salt flux at each level (dₜ∫SdV)"))
        defVar(new_ds, "∫dₜ∫SdVdV", Float64, tuple("time_derivative"),
                attrib = Dict("longname" => "Volume integrated salt flux (∫dₜ∫SdVdV)"))
        defVar(new_ds, "dₜ∫Sz✶dV", Float64, tuple("time_derivative"),
                attrib = Dict("longname" => "Volume integrated background salinity potential energy (dₜ∫Sz✶dV)."))

        for t ∈ eachindex(time)
            S = reshape(ds[:S][:, :, :, t], :)
            p = sortperm(S)
            new_ds[:∫Sz✶dV][t] = sum(S .* z[p] * ΔV)
            new_ds[:∫SdV][:, t] = cumsum(S[p] * ΔV)
        end

        for t ∈ 1:length(time)-1
            new_ds[:dₜ∫SdV][:, t] = vec(diff(new_ds[:∫SdV][:, t:t+1], dims = 2)) / Δt[t]
            new_ds[:∫dₜ∫SdVdV][t] = sum(new_ds[:dₜ∫SdV][:, t] * ΔV)
            new_ds[:dₜ∫Sz✶dV][t] = diff(new_ds[:dₜ∫SdV][t:t+1])[1] / Δt[t]
        end

    end

    close(ds)

    return nothing
end

tracers = "tracers.nc"
file = "salt_checks.nc"
salt_checks!(file, tracers)
