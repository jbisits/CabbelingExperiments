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

    Δt = diff(time)
    Δx = diff(ds[:xC][1:2])[1]
    Δy = diff(ds[:yC][1:2])[1]
    Δz = diff(ds[:zC][1:2])[1]

    ΔV = Δx * Δy * abs(Δz)

    # Sorting
    Csorted = Array{Float64}(undef, length(reshape(C[:, :, :, 1], :)), length(time))
    ∫CdV = similar(Csorted)
    ∂z_Csorted = Array{Float64}(undef, length(reshape(C[:, :, :, 1], :))-1, length(time))

    for t ∈ eachindex(time)

        Csorted[:, t] = sort(reshape(C[:, :, :, t], :))
        ∫CdV[:, t] = cumsum(Csorted[:, t] * ΔV)
        ∂z_Csorted[:, t] = diff(Csorted[:, t]) / Δz

    end

    ∂z_Csorted_mean = mean(∂z_Csorted, dims = 2)

    # Tracer flux
    dₜ∫CdV = Array{Float64}(undef, length(reshape(C[:, :, :, 1], :)), length(time)-1)

    for t ∈ 1:length(time)-1

        dₜ∫CdV[:, t] = vec(diff(∫CdV[:, t:t+1], dims = 2)) / Δt[t]

    end

    κ_forward = Array{Float64}(undef, length(reshape(C[:, :, :, 1], :))-1, length(time)-1)
    κ_mean = Array{Float64}(undef, length(reshape(C[:, :, :, 1], :))-1, length(time)-1)

    κ_forward = dₜ∫CdV[2:end, :] ./ ∂z_Csorted[:, 2:end]
    for t ∈ 1:length(time)-1

        κ_mean[:, t] = dₜ∫CdV[2:end, t] ./ ∂z_Csorted_mean

    end

    V = eachindex(Csorted[:, 1]) * ΔV
    SA = 0.1 * 0.1 # Surface area, would be better if not hard coded
    equivalent_z = -V / SA

    close(ds)

    NCDataset(filename, "a") do ds

        defDim(ds, "cumulative_volume", length(V))
        defDim(ds, "cumulative_volume_derivative", length(V)-1)
        defDim(ds, "time", length(time))
        defDim(ds, "time_derivative", length(Δt))
        defVar(ds, "volume", V, tuple("cumulative_volume"))
        defVar(ds, "equivalent_z", equivalent_z, tuple("cumulative_volume"))
        defVar(ds, "time", time, tuple("time"))
        defVar(ds, "time_derivative", cumsum(Δt), tuple("time_derivative"))
        defVar(ds, string.(tracer) * "sorted", Csorted, ("cumulative_volume", "time"),
               attrib = Dict("longname" => "Sorted (smallest to largest) tracer field " * string.(tracer)))
        defVar(ds, "∂z_" * string.(tracer) * "sorted", ∂z_Csorted, ("cumulative_volume_derivative", "time"),
               attrib = Dict("longname" => "Vertical derivative of sorted tracer field " * string.(tracer)))
        defVar(ds, "∂z_" * string.(tracer) * "sorted_mean", ∂z_Csorted_mean, tuple("cumulative_volume_derivative"),
               attrib = Dict("longname" => "Time mean of vertical derivative of sorted tracer field " * string.(tracer)))
        defVar(ds, "dₜ∫" * string.(tracer) * "dV", dₜ∫CdV, tuple("cumulative_volume", "time_derivative"),
               attrib = Dict("longname" => "Vertical flux of " * string.(tracer) * " into a fixed volume (i.e. a level)"))
        defVar(ds, "κ_effective" * string.(tracer), κ_forward, tuple("cumulative_volume_derivative", "time_derivative"),
               attrib = Dict("longname" => "Effective vertical diffusivity for " * string.(tracer)))
        defVar(ds, "κ_effective_timemean" * string.(tracer), κ_mean, tuple("cumulative_volume_derivative", "time_derivative"),
               attrib = Dict("longname" => "Effective vertical diffusivity for " * string.(tracer) * "using a time mean of the vertical derivative"))


    end

    return nothing
end

tracers = "tracers.nc"
compute_diffusivity(tracers, :T)
compute_diffusivity(tracers, :S)
