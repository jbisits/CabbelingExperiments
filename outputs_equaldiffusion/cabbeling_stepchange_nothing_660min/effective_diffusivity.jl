# Script to compute an effective vertical diffusivity from a `sort`ed tracer field using the
# tracer percentile framework (Sohail et al. (2021)).
using NCDatasets, StatsBase

"""
    function makefile(filename::AbstractString, saved_output::AbstractString)
Create a file with `filename` that has the same dimensions as `saved_output`.
"""
function makefile(filename::AbstractString, saved_output::AbstractString)

    co = NCDataset(saved_output)

    NCDataset(joinpath(@__DIR__, filename), "c") do ds

        for key ∈ keys(co.attrib)
            ds.attrib[key] = co.attrib[key]
        end

    end
    close(co)

    return nothing

end
"""
    function compute_diffusivity(saved_output, tracer)
Compute the effective vertical diffusivity for `tracer` in `saved_output`. The effective
difffusivity is calculated using the temperature percentile framework (Sohail et al. (2021)).
First the tracer content within a fixed volume is calculated by `sort`ing the tracer values.
The tracer content tendency is then the time_derivative of the (fixed) volume integrated
tracer content. The effective diffusivity is found by dividing the tendency by the vertical
tracer gradient.
"""
function compute_diffusivity(saved_output::AbstractString, tracer::Symbol)

    filename = string.(tracer) * "_effective_diffusivity.nc"
    isfile(filename) ? nothing : makefile(filename, saved_output)

    ds = NCDataset(saved_output)

    time = ds[:time][:]
    C = ds[tracer]
    C_size = size(C)
    C_length = prod(C_size[1:3])

    Δt = diff(time)
    Δx = diff(ds[:xC][1:2])[1]
    Δy = diff(ds[:yC][1:2])[1]
    Δz = diff(ds[:zC][1:2])[1]

    ΔV = Δx * Δy * abs(Δz)

    V = (1:C_length) * ΔV
    SA = 0.1 * 0.1 # Surface area, would be better if not hard coded
    equivalent_z = -V / SA

    NCDataset(filename, "a") do ds_diff

        defDim(ds_diff, "cumulative_volume", length(V))
        defDim(ds_diff, "cumulative_volume_derivative", length(V)-1)
        defDim(ds_diff, "time", length(time))
        defDim(ds_diff, "time_derivative", length(Δt))
        defVar(ds_diff, "volume", V, tuple("cumulative_volume"))
        defVar(ds_diff, "equivalent_z", equivalent_z, tuple("cumulative_volume"))
        defVar(ds_diff, "time", time, tuple("time"))
        defVar(ds_diff, "time_derivative", cumsum(Δt), tuple("time_derivative"))
        defVar(ds_diff, string.(tracer) * "sorted", Float64, ("cumulative_volume", "time"),
               attrib = Dict("longname" => "Sorted (smallest to largest) tracer field " * string.(tracer)))
        defVar(ds_diff, "∫" * string.(tracer) * "dV", Float64, ("cumulative_volume", "time"),
               attrib = Dict("longname" => "Vertical flux for " * string.(tracer)))
        defVar(ds_diff, "∂z_" * string.(tracer) * "sorted", Float64, ("cumulative_volume_derivative", "time"),
               attrib = Dict("longname" => "Vertical derivative of sorted tracer field " * string.(tracer)))

        for t ∈ eachindex(time)

            C_sorted = sort(reshape(C[:, :, :, t], :))
            ds_diff[string.(tracer) * "sorted"][:, t] = C_sorted
            ds_diff["∫" * string.(tracer) * "dV"][:, t] = cumsum(C_sorted * ΔV)
            ds_diff["∂z_" * string.(tracer) * "sorted"][:, t] = diff(C_sorted) / Δz

        end

    end

    close(ds)

    NCDataset(filename, "a") do ds_diff

        defVar(ds_diff, "dₜ∫" * string.(tracer) * "dV", Float64, tuple("cumulative_volume", "time_derivative"),
               attrib = Dict("longname" => "Vertical flux of " * string.(tracer) * " into a fixed volume (i.e. a level)"))
        defVar(ds_diff, "κ_effective" * string.(tracer), Float64, tuple("cumulative_volume_derivative", "time_derivative"),
               attrib = Dict("longname" => "Effective vertical diffusivity for " * string.(tracer)))

        for t ∈ 1:length(time)-1

            flux = vec(diff(ds_diff["∫" * string.(tracer) * "dV"][:, t:t+1], dims = 2)) / Δt[t]
            ds_diff["dₜ∫" * string.(tracer) * "dV"][:, t] = flux
            ds_diff["κ_effective" * string.(tracer)][:, t] = flux[2:end] ./ ds_diff["∂z_" * string.(tracer) * "sorted"][:, t+1]

        end

    end

    return nothing
end

tracers = "tracers.nc"
compute_diffusivity(tracers, :T)
compute_diffusivity(tracers, :S)
