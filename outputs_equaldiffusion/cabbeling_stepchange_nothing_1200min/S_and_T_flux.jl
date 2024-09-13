using NCDatasets, JLD2, GibbsSeaWater

tracers = "tracers.nc"
velocities = "velocities.nc"
bflux = "buoyancy_flux_interp_face.jld2"

ds_tracers = NCDataset(tracers)
time = ds_tracers[:time][:]
ΔV = diff(ds_tracers[:xC][1:2])[1] * diff(ds_tracers[:yC][1:2])[1] * diff(ds_tracers[:zC][1:2])[1]
ds_vel = NCDataset(velocities)

∫Sw = similar(time)
∫Θw = similar(time)
for t ∈ eachindex(time)

    S = ds_tracers[:S][:, :, :, t]
    S1 = @view S[:, :, 1:end-1]
    S2 = @view S[:, :, 2:end]
    S_interp = cat(S[:, :, 1], 0.5 * (S1 .+ S2), S[:, :, end], dims = 3)
    T = ds_tracers[:T][:, :, :, t]
    T1 = @view T[:, :, 1:end-1]
    T2 = @view T[:, :, 2:end]
    T_interp = cat(T[:, :, 1], 0.5 * (T1 .+ T2), T[:, :, end], dims = 3)
    w = ds_vel[:w][:, :, :, t]

    ∫Sw[t] = sum(S_interp .* w) * ΔV
    ∫Θw[t] = sum(T_interp .* w) * ΔV

end

close(ds_tracers)
close(ds_vel)

jldopen(bflux, "a+") do file
    file["∫Sw"] = ∫Sw
    file["∫Θw"] = ∫Θw
    file["g"] = 9.81
end
