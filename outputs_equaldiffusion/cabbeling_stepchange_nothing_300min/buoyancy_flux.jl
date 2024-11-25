using NCDatasets, JLD2

co = "computed_output.nc"
velocities = "velocities.nc"
bflux = "buoyancy_flux_interp_face.jld2"

ds_co = NCDataset(co)
time = ds_co[:time][:]
ΔV = diff(ds_co[:xC][1:2])[1] * diff(ds_co[:yC][1:2])[1] * diff(ds_co[:zC][1:2])[1]
ds_vel = NCDataset(velocities)

g = 9.81
∫gρw = similar(time)
for t ∈ eachindex(time)

    σ = ds_co[:σ][:, :, :, t]
    σ1 = @view σ[:, :, 1:end-1]
    σ2 = @view σ[:, :, 2:end]
    σ_interp = cat(σ[:, :, 1], 0.5 * (σ1 .+ σ2), σ[:, :, end], dims = 3)
    w = ds_vel[:w][:, :, :, t]

    ∫gρw[t] = g * sum(σ_interp .* w) * ΔV

end
close(ds_co)
close(ds_vel)

NCDataset(co, "a") do ds_co
    defVar(ds_co, "∫gρw", ∫gρw, ("time",),
            attrib = Dict("longname" => "Volume integrated density flux in post processing"))
end

jldopen(bflux, "w") do file
    file["time"] = time
    file["∫gρw"] = ∫gρw
end
