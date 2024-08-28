using NCDatasets, JLD2

co = "computed_output.nc"
velocities = "velocities.nc"
bflux = "buoyancy_flux.jld2"

ds_co = NCDataset(co)
time = ds_co[:time][:]
ΔV = diff(ds_co[:xC][1:2])[1] * diff(ds_co[:yC][1:2])[1] * diff(ds_co[:zC][1:2])[1]
ds_vel = NCDataset(velocities)

g = 9.81
∫gρw = similar(time)
for t ∈ eachindex(time)

    σ = ds_co[:σ][:, :, :, t]
    w = ds_vel[:w][:, :, :, t]

    ∫gρw[t] = g * sum(σ .* w) * ΔV

end

close(ds_co)
close(ds_vel)

jldopen(bflux, "w") do file
    file["time"] = time
    file["∫gρw"] = ∫gρw
end
