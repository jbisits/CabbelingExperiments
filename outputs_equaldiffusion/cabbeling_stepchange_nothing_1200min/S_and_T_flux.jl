using NCDatasets, JLD2, GibbsSeaWater

tracers = "tracers.nc"
velocities = "velocities.nc"
bflux = "buoyancy_flux.jld2"

ds_tracers = NCDataset(tracers)
time = ds_tracers[:time][:]
ΔV = diff(ds_tracers[:xC][1:2])[1] * diff(ds_tracers[:yC][1:2])[1] * diff(ds_tracers[:zC][1:2])[1]
ds_vel = NCDataset(velocities)

∫βSw = similar(time)
∫αΘw = similar(time)
for t ∈ eachindex(time)

    S = ds_tracers[:S][:, :, :, t]
    T = ds_tracers[:T][:, :, :, t]
    α = gsw_alpha.(S, T, 0)
    β = gsw_beta.(S, T, 0)
    w = ds_vel[:w][:, :, :, t]
    w1 = @view w[:, :, 1:end-1]
    w2 = @view w[:, :, 2:end]
    w_interp = (w1 .+ w2) / 2

    ∫βSw[t] = sum(β .* S .* w_interp) * ΔV
    ∫αΘw[t] = sum(α .* T .* w_interp) * ΔV

end

close(ds_tracers)
close(ds_vel)

jldopen(bflux, "a+") do file
    file["∫βSw"] = ∫βSw
    file["∫αΘw"] = ∫αΘw
    file["g"] = 9.81
end
