using NCDatasets, GibbsSeaWater

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

    copy_vars = ("volume", "equivalent_z", "time", "time_derivative")
    NCDataset(filename, "a") do ds

        for v ∈ copy_vars
            defVar(ds, co[v])
        end

    end

    close(co)

    return nothing

end

"""
    function compute_buoyancy_flux(bflux_file, T_output, S_output)
Compute the bouyancy flux from the temperature and salinity fluxes saved in `T_output` and
`S_output` and save (with some other variables) in `bflux_file`.
"""
function compute_buoyancy_flux!(bflux_file::AbstractString, T_output::AbstractString, S_output::AbstractString)

    ds_S, ds_T = NCDataset(S_output), NCDataset(T_output)
    time_derivative = ds_T[:time_derivative][:]
    ΔV = diff(ds_T[:cumululative_volume][1:2])[1]
    cₚ = 4000 # specific heat capacity
    g = 9.81
    ds_bflux = NCDatasets(bflux_file, "a")
    defVar(ds_bflux, "J_b", Float64, tuple("cumulative_volume", "time_derivative"),
            attrib = Dict("longname" => "Bouyancy flux, J = -g(αFₜ/Cₚ - βFₛ), computed from T and S fluxes."))
    defVar(ds_∫bflux, "∫J_b", Float64, tuple("time_derivative"),
            attrib = Dict("longname" => "Volume integrated bouyancy flux, ∫JdV = -g∫(αFₜ/Cₚ - βFₛ)dV, computed from T and S fluxes."))

    for t ∈ eachindex(time_derivative)

        S_sorted, T_sorted = ds_S[:S_sorted][:], ds_T[:T_sorted][:]
        α, β = gsw_alpha.(S_sorted, T_sorted, 0), gsw_beta.(S_sorted, T_sorted, 0)
        Fₛ, Fₜ = ds_S[:dₜ∫SdV], ds_T[:dₜ∫TdV]
        ds_bflux[:, t] = @. -g * ((α * Fₜ) / cₚ - β * Fₛ)
        ds_∫bflux[t] = sum(ds_bflux[:, t]) * ΔV

    end

    return nothing

end

T_output_diff = "T_effective_diffusivity.nc"
S_output_diff = "S_effective_diffusivity.nc"
bflux_file = "buoyancy_flux.nc"
isfile(bflux_file) ? nothing : makefile(bflux_file, T_output_diff)
compute_buoyancy_flux!(bflux_file, T_output_diff, S_output_diff)
