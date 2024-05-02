### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a290de4c-f791-11ee-0e04-11e4e34b7428
begin
	using Pkg
	Pkg.activate("..")
	using JLD2, Oceananigans, CairoMakie, PlutoUI, Statistics
end

# ╔═╡ ac639feb-9ee4-43f4-acd9-466fe3478d40
begin
	choose_expt = @bind experiment Select(["isothermal_nonoise", "isothermal_withnoise", "cabbeling_cd1_nonoise", "cabbeling_cd1_withnoise", "cabbeling_cd10_nonoise", "cabbeling_cd10_withnoise"])
	
	md"""
	# Flux and background state check with 1D model
	
	While computing the background state and the buoyancy flux for the model there were some discrepancies.
	In this notebook I run equivalent 1D model and analysis to try and find out what might be going on with some of these calculations.
	
	The main two experiments I will compare are the isothermal experiment and the cabbeling experiment (this is all I have done for DNS).
	
	## Model runs
	
	I have three model runs with background diffusivity ``\kappa_{back} = 1\times 10^{-2}\mathrm{m}^{2}\mathrm{s}^{-1}``:
	
	1. isothermal, uniform ``T = 0.5^{\circ}C`` (there is also a lower resolution one of these experiments)
	2. cabbeling with convective diffusivity ``\kappa_{c} = 1\mathrm{m}^{2}\mathrm{s}^{-1}``
	2. cabbeling with convective diffusivity ``\kappa_{c} = 10\mathrm{m}^{2}\mathrm{s}^{-1}``
	
	In all experiments with no initial noise in the salinity field, the salinity profile and the background profile *should be the same* as there are no sources and sinks of salinity.
	We see this in the output below.
	
	In the isothermal case, the background state (for density) should be equal to the profile so the BPE and the PE should be equal (or very close to equal) throughout the length of the because the mixed water at the interface will not be denser/lighter than the lower/upper layers.
	This means that there is no available potential energy in the system.
	
	In the cabbeling case, the background state (for density) and profile should start roughly equal but after the mixing has created some denser water the PE and the BPE should separate until when/if lower layer has all been transformed to the maximum density.
	This means that there is available potential energy in the system which is dissipated as the simulation runs.

	Adding noise about the interface gives a closer representation of what the DNS experiments look like.

	The different experiments can be chose from this list $(choose_expt).
	"""
end

# ╔═╡ edd033e9-55ca-4b68-bf03-330baaa35e67
md"""
## Salinity
"""

# ╔═╡ d2a1ec90-16ef-47b2-9adb-8829681ad5d9
begin
	iso_data = "../1DModel/OneDModelOutput_"*experiment*".jld2"
	S = FieldTimeSeries(iso_data, "S")
	T = FieldTimeSeries(iso_data, "T")
	σ₀ = FieldTimeSeries(iso_data, "σ")
	t = S.times
	Δt = diff(t)
	zplot = znodes(S)
	z = reverse(abs.(znodes(S)))
	Δz = S.grid.Δzᵃᵃᶜ

	κc = if experiment ∈ ("cabbeling_cd1_nonoise", "cabbeling_cd1_withnoise")
			 1.0
		 else
			 10.0
		 end

	#Computations to check

	# Potential energy for salt
	∫Szdz = sum(interior(S, 1, 1, :, :) .* z * Δz, dims = 1)
	dₜ∫Szdz = vec(diff(∫Szdz, dims = 2)) ./ Δt

	# Salt sorted, content and flux
	S✶ = similar(S.data[1, 1, :, :])
	∫Sdz = similar(S.data[1, 1, :, :])
	for i ∈ eachindex(t)
		S✶[:, i] = sort(interior(S, 1, 1, :, i), rev = true)
		∫Sdz[:, i] = cumsum(reverse(S✶[:, i]) * Δz)
	end

	dₜ∫Sdz = diff(∫Sdz, dims = 2)./ Δt'
	∫dₜ∫Sdzdz = vec(sum(dₜ∫Sdz * Δz, dims = 1))

	# Background potential energy
	∫S✶zdz = sum(S✶ .* z * Δz, dims = 1)
	dₜ∫S✶zdz = vec(diff(∫S✶zdz, dims = 2)) ./ Δt

	nothing
end

# ╔═╡ 0efe7042-82d9-44aa-be55-eb3634425bd5
let
	fig, ax, hm = heatmap(t, zplot, interior(S, 1, 1, :, :)', colormap = :haline)
	ax.title = "Salinity hovmoller"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm, label = "S (g/kg )")
	fig
end

# ╔═╡ daf8bad7-444b-4da7-968b-32f29f1b7eba
md"""
### Time evolution
"""

# ╔═╡ 06bd6b0a-8d01-43c1-8560-f3c20972fe9e
begin
	t_slider = @bind timestep PlutoUI.Slider(eachindex(t))
	md"""
	Slider for plotting salinity profiles (profile and sorted profile)
	$(t_slider)
	"""
end

# ╔═╡ 08a90e73-8312-4cd2-9c5d-8ff9829962bc
let
	fig, ax = lines(interior(S, 1, 1, :, timestep), zplot, label = "Profile")
	lines!(ax, S✶[:, timestep], zplot, label = "Sorted profile", linestyle = :dash)
	ax.title = "Salinity profiles t = $(t[timestep] / 60) minutes"
	ax.xlabel = "Salinity (g/kg)"
	ax.ylabel = "z (m)"
	axislegend(ax, position = :lb)
	fig
end

# ╔═╡ ebff4d86-8013-4473-ae43-0fc2b4fe2015
md"""
### "Energetics" (salinity only)

Plot of potential, background potential and available potential energies using only salt tracer (so not really energetic quantities...).
The quantities (without gravity) are calculated as
```math
\begin{aligned}
\frac{\mathrm{d}}{\mathrm{d}t}E_{p} &= \frac{\mathrm{d}}{\mathrm{d}t}\int S z dz \\
\frac{\mathrm{d}}{\mathrm{d}t}E_{b} &= \frac{\mathrm{d}}{\mathrm{d}t}\int S^{*} z dz \\
\frac{\mathrm{d}}{\mathrm{d}t}E_{a} &= \frac{\mathrm{d}}{\mathrm{d}t}E_{p} - \frac{\mathrm{d}}{\mathrm{d}t}E_{b}
\end{aligned}
```
where ``S^{*}`` is the sorted salinity profile.
"""

# ╔═╡ 1b7fba62-557e-4c53-8f67-b45b16548e1f
let
	fig2 = Figure(size = (500, 1000))
	ax1 = Axis(fig2[1, 1], title = "Salinity potential energies")
	ylims!(ax1, maximum(∫Szdz) .+ [-1e4, 1e3])
	lines!(ax1, t, vec(∫Szdz), label = "PE")
	lines!(ax1, t, vec(∫S✶zdz), label = "BPE", linestyle = :dash)
	axislegend(ax1, position = :rb)
	ax2 = Axis(fig2[2, 1], title = "Available potential energy")
	lines!(ax2, t, vec(∫Szdz) .- vec(∫S✶zdz), label = "APE")
	axislegend(ax2, position = :rt)
	fig2
end

# ╔═╡ 944526a6-5154-44ac-af06-c7dba206007c
begin
	fig2 = Figure(size = (500, 1000))
	ax1 = Axis(fig2[1, 1], title = "Time change salinity potential energies")
	lines!(ax1, t[2:end], dₜ∫Szdz, label = "dₜPE")
	lines!(ax1, t[2:end], dₜ∫S✶zdz, label = "dₜBPE", linestyle = :dash)
	axislegend(ax1, position = :rb)
	ax2 = Axis(fig2[2, 1], title = "Time change available potential energy")
	lines!(ax2, t[2:end], dₜ∫Szdz .- dₜ∫S✶zdz, label = "dₜAPE", color = :red)
	axislegend(ax2, position = :rt)
	fig2
end

# ╔═╡ 3e3eb9fd-8208-4f99-8cbb-730e4a501772
md"""
### Flux and BPE

Changes to the background potential energy can **only occur due to vertical diffusive flux** so we should have equality between
```math
\frac{\mathrm{d}}{\mathrm{d}t}E_{b} \quad \mathrm{and} \quad \frac{\mathrm{d}}{\mathrm{d}t}\iint S^{*} dz dz
```
"""

# ╔═╡ ab6e5cd8-9a90-4cc2-8e83-330aa7ae8b30
let
	fig = Figure(size = (500, 1000))
	ax = Axis(fig[1, 1], title = "Salt flux and BPE")
	lines!(ax, t[2:end], dₜ∫S✶zdz, label = "dₜ∫S✶zdz")
	lines!(ax, t[2:end], ∫dₜ∫Sdzdz, label = "∫dₜ∫Sdzdz", linestyle = :dash)
	axislegend(ax, position = :rb)
	ax2 = Axis(fig[2, 1],
				title = "Absolute error between salinity flux and background salt energy",
				xlabel = "time (s)",
				ylabel = "Absolute error")
	lines!(ax2, t[2:end], abs.(∫dₜ∫Sdzdz - dₜ∫S✶zdz))
	fig
end

# ╔═╡ 16562398-cd61-4697-b94a-931a5d44a2f3
md"""
## Density

As the salinity is looking ok I would like to look at density and see if I can make sense of things in the 1D case.
"""

# ╔═╡ b3fc213d-03c1-4fb5-9c13-f3f6e5f1f648
begin
	# Potential energy no gravity
	g = 1 # 9.81
	∫σzdz = g * sum(interior(σ₀, 1, 1, :, :) .* z * Δz, dims = 1)
	dₜ∫σzdz = vec(diff(∫σzdz, dims = 2)) ./ Δt
 
	# Sorted density, density contet, density flux
	σ✶ = similar(σ₀.data[1, 1, :, :])
	∫σdz = similar(σ₀.data[1, 1, :, :])
	for i ∈ eachindex(t)
		σ✶[:, i] = sort(interior(σ₀, 1, 1, :, i), rev = true)
		∫σdz[:, i] = cumsum(reverse(σ✶[:, i]) * Δz)
	end

	dₜ∫σdz = diff(∫σdz, dims = 2)
	∫dₜ∫σdzdz = vec(sum(dₜ∫σdz * Δz, dims = 1)) ./ Δt

	# Background potential energy
	∫σ✶zdz = g * sum(σ✶ .* z * Δz, dims = 1)
	dₜ∫σ✶zdz = vec(diff(∫σ✶zdz, dims = 2)) ./ Δt

	nothing
end

# ╔═╡ 385c06a8-f60b-4f5b-8bd2-e2a56447397c
md"""
### Time evolution
"""

# ╔═╡ 2a164d51-759d-4341-811d-d1fd406c0c3d
let
	fig, ax, hm = heatmap(t, zplot, interior(σ₀, 1, 1, :, :)', colormap = :dense)
	ax.title = "Density hovmoller"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm, label = "σ₀ (kg/m^3 )")
	fig
end

# ╔═╡ 6cce82f4-b6a8-4be4-8e15-e93908d72eee
begin
	t_slider2 = @bind timestep2 PlutoUI.Slider(eachindex(t))
	md"""
	Slider for plotting density profiles (profile and sorted profile)
	$(t_slider2)
	"""
end

# ╔═╡ febc3ebe-8d07-4c7d-aaec-fd04bb7378a2
let
	fig, ax = lines(interior(σ₀, 1, 1, :, timestep2), zplot, label = "Initial profile")
	lines!(ax, σ✶[:, timestep2], zplot, label = "Sorted profile", linestyle = :dash)
	ax.title = "Density profiles t = $(t[timestep2] / 60) minutes"
	ax.xlabel = "Density (σ₀, kg/m^3)"
	ax.ylabel = "z (m)"
	axislegend(ax, position = :lb)
	fig
end

# ╔═╡ c98ca6ac-0af9-4bbf-88d7-8ef72a58b3e0
md"""
### Energetics
"""

# ╔═╡ 4fea18d5-9fb5-40ff-8ea0-7e9bbf3cb9e1
let
	fig2 = Figure(size = (500, 1000))
	ax1 = Axis(fig2[1, 1], title = "Potential energies")
	ylims!(ax1, maximum(∫σzdz) .+ [-1e2, 1e2])
	lines!(ax1, t, vec(∫σzdz), label = "PE")
	lines!(ax1, t, vec(∫σ✶zdz), label = "BPE", linestyle = :dash)
	axislegend(ax1, position = :rb)
	ax2 = Axis(fig2[2, 1], title = "Available potential energy")
	lines!(ax2, t, vec(∫σzdz) .- vec(∫σ✶zdz), label = "APE")
	axislegend(ax2, position = :rt)
	fig2
end

# ╔═╡ 3217dc82-a988-448f-9bd4-ef5f62b75630
let
	fig2 = Figure(size = (500, 1000))
	ax1 = Axis(fig2[1, 1], title = "Time change potential energies")
	lines!(ax1, t[2:end], dₜ∫σzdz, label = "dₜPE")
	lines!(ax1, t[2:end], dₜ∫σ✶zdz, label = "dₜBPE", linestyle = :dash)
	axislegend(ax1, position = :rb)
	ax2 = Axis(fig2[2, 1], title = "Time change available potential energy")
	lines!(ax2, t[2:end], dₜ∫σzdz .- dₜ∫σ✶zdz, label = "dₜAPE", color = :red)
	axislegend(ax2, position = :rt)
	fig2
end

# ╔═╡ e01432d2-627c-44cf-aba7-67d2094092d3
md"""
### Flux and BPE
"""

# ╔═╡ 1a9eee8b-545f-4da3-a931-49d447521e3c
let
	fig = Figure(size = (500, 1000))
	ax = Axis(fig[1, 1], title = "Density flux and BPE")
	lines!(ax, t[2:end], dₜ∫σ✶zdz, label = "dₜ∫σ✶zdz")
	lines!(ax, t[2:end], ∫dₜ∫σdzdz, label = "∫dₜ∫Sdzdz", linestyle = :dash)
	axislegend(ax, position = :rb)
	ax2 = Axis(fig[2, 1],
				title = "Absolute error between density flux and background potential energy",
				xlabel = "time (s)",
				ylabel = "Absolute error (log10)")
	lines!(ax2, t[2:end], log10.(abs.(∫dₜ∫σdzdz - dₜ∫σ✶zdz)))
	fig
end

# ╔═╡ 0121a899-a1bf-4fc8-9aa2-ce5c801753ec
md"""
## Diffusivity

### Salt (or temperature)

We calculate the diffusivity from the tracers using
```math
\kappa_{S} = \frac{\mathrm{d}_{t}\int S^{*} \mathrm{d}z}{\mathrm{d}S^{*} / \mathrm{d}z}
```

As this model has the diffusivity values parameterised we should recover the background diffusivity ``\kappa_{back} = 1\times 10^{-2}m^{2}s^{-1}`` in the isothermal case and something that looks like the convective adjustment schemes used above.

From output below can see for all cases there diffusivity about the estimates match what is expected --- background diffusivity (``\kappa_{back} = 1\times 10^{-2}m^{2}s^{-1}``) in isothermal case convective diffusivity ``(\kappa_{c} = 1m^{2}s^{-1}``, ``\kappa_{c} = 10m^{2}s^{-1})`` in cabbeling cases.
"""

# ╔═╡ 9a82d299-0274-4e68-9c4b-da1350e52fe1
begin
	dSdz = diff(reverse(S✶, dims = 1), dims = 1) ./ Δz
	replace!(dSdz, 0 => NaN)
	κₛ = dₜ∫Sdz[2:end, :] ./ dSdz[:, 2:end]
	zrange = 690:710
	∫κₛ = sum(κₛ[.!isnan.(κₛ)] * Δz, dims = 1)
	κₛ_mean = mean(κₛ[.!isnan.(κₛ)])
	nothing
end

# ╔═╡ 2642bc19-9c7c-45f5-b199-fe318da98101
begin
	z_slider = @bind z_ PlutoUI.Slider(eachindex(z[1:end-1]))
	md"""
	Choose a depth level to display the diffusivity calculated $(z_slider).
	"""
end

# ╔═╡ 52a15ea8-600d-4bbf-8dee-def4895a4ded
let
	# zrange = 680:720
	# fig, ax, hm = heatmap(t[2:end], z[zrange], κₛ[zrange, :]')
	# Colorbar(fig[1, 2], hm)
	# fig
	# fig, ax = series(t[20:end-1], κₛ[zrange, 20:end], labels = ["z = -$(round(i, digits = 2))m" for i ∈ z[zrange]], solid_color=:black)
	# ax.title = "Diffusivity estimates for salinity about the model interface"
	# ax.xlabel = "time (s)"
	# ax.ylabel = "κ (m2/s)"
	# axislegend(ax)
	# Legend(fig[1, 2], ax)
	fig, ax = lines(t[20:end-1], κₛ[z_, 20:end], label = "z = -$(round(z[z_], digits = 2))", solid_color=:black)
	ax.title = "Diffusivity estimates for salinity about the model interface"
	ax.xlabel = "time (s)"
	ax.ylabel = "κ (m2/s)"
	hlines!(0.01, label = "κ background", color = :black, linestyle = :dash)
	hlines!(κc, label = "κ convective", color = :red, linestyle = :dash)
	axislegend(ax)
	fig
end

# ╔═╡ cd11a001-c306-4e16-8df5-8a51eb0a46ee
let
	zrange = 690:710
	# fig, ax, hm = heatmap(t[2:end], z[zrange], κₛ[zrange, :]')
	# Colorbar(fig[1, 2], hm)
	# fig
	fig, ax = series(t[20:end-1], κₛ[zrange, 20:end], labels = ["z = -$(round(i, digits = 2))m" for i ∈ z[zrange]], solid_color=:grey)
	ax.title = "Diffusivity estimates for salinity about the model interface"
	ax.xlabel = "time (s)"
	# ax.ylabel = "κ (m2/s)"
	# axislegend(ax)
	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 4f5f7717-53fc-4252-804a-c0a3293df8d5
begin
	iso_path = "../outputs_equaldiffusion/isothermal_stepchange_nothing_660min/isothermal_profile.jld2"
	cab_path = "../outputs_equaldiffusion/cabbeling_stepchange_nothing_660min/cabbeling_profile.jld2"
	DNS_EXPT = @bind dns_expt Select([iso_path => "isothermal", cab_path => "cabbeling"])
	md"""
	# DNS profile
	
	I am now interested in how a single profile of the DNS compares to this.
	Eventually I will also look at a horizontally averaged profile but the full thing will use all volume elements of the 3D grid resorted.
	
	There are only two DNS experiments one isothermal, single scalar field (salinity), and the cabbeling, two tracer fields (salinity and temperature).
	
	There are the same number of vertical levels so I should be able as the 1D model above but the diffusivity is a `ScalarDiffusivity` with value of ``1\times 10^{-7}\mathrm{m}^{2}\mathrm{s}^{-1}``.

	Simulation: $(DNS_EXPT).

	## Salinity
	"""
end

# ╔═╡ 0602c590-4435-4963-acd3-15798c9de894
begin
	file = jldopen(dns_expt)
	S_dns = file["S"]
	T_dns = file["T"]
	σ₀_dns = file["σ"]
	t_dns = file["time"]
	close(file)
	
	Δt_dns = diff(t_dns)

	Δz_dns = dns_expt == cab_path ? 0.0007142857142857784 :  0.0010000000000000009
	z_dns_plot = range(-1, 0, step = Δz_dns)
	z_dns = reverse(abs.(z_dns_plot))
	
	#Computations to check

	# Potential energy for salt
	∫Szdz_dns = sum(S_dns .* z_dns * Δz_dns, dims = 1)
	dₜ∫Szdz_dns = vec(diff(∫Szdz_dns, dims = 2)) ./ Δt_dns

	# Salt sorted, content and flux
	S✶_dns = similar(S_dns)
	∫Sdz_dns = similar(S_dns)
	for i ∈ eachindex(t_dns)
		S✶_dns[:, i] = sort(S_dns[:, i], rev = true)
		∫Sdz_dns[:, i] = cumsum(reverse(S✶_dns[:, i]) * Δz_dns)
	end

	dₜ∫Sdz_dns = diff(∫Sdz_dns, dims = 2)./ Δt_dns'
	∫dₜ∫Sdzdz_dns = vec(sum(dₜ∫Sdz_dns * Δz_dns, dims = 1))

	# Background potential energy
	∫S✶zdz_dns = sum(S✶_dns .* z_dns * Δz_dns, dims = 1)
	dₜ∫S✶zdz_dns = vec(diff(∫S✶zdz_dns, dims = 2)) ./ Δt_dns

	nothing
end

# ╔═╡ ebf2aa23-5cdd-4de7-a285-ce1534666899
md"""
### Time evolution
"""

# ╔═╡ 7a2b5232-e081-45ae-95be-7d219aae91ed
let
	fig, ax, hm = heatmap(t_dns, z_dns, S_dns', colormap = :haline)
	ax.title = "Salinity hovmoller"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm, label = "S (g/kg )")
	fig
end

# ╔═╡ ffdd23ec-190f-4307-adc9-7058bfcc2593
begin
	t_slider_dns = @bind timestep_dns PlutoUI.Slider(eachindex(t_dns))
	md"""
	Slider for plotting salinity profiles (profile and sorted profile)
	$(t_slider_dns)
	"""
end

# ╔═╡ b4dcdf1e-09ce-41e5-83de-90356dce3eeb
let
	fig, ax = lines(S_dns[:, timestep_dns], z_dns_plot, label = "Profile")
	lines!(ax, S✶_dns[:, timestep_dns], z_dns_plot, label = "Sorted profile", linestyle = :dash)
	ax.title = "Salinity profiles t = $(t[timestep_dns] / 60) minutes"
	ax.xlabel = "Salinity (g/kg)"
	ax.ylabel = "z (m)"
	axislegend(ax, position = :lb)
	fig
end

# ╔═╡ 7579eacb-fddd-40bb-a50f-a8f98cbbdbcd
md"""
### "Energetics" (Salinity only) 
"""

# ╔═╡ cf9aba94-ef48-4478-a49a-92e9da84e9d9
let
	fig2 = Figure(size = (500, 1000))
	ax1 = Axis(fig2[1, 1], title = "Salinity potential energies")
	lines!(ax1, t_dns, vec(∫Szdz_dns), label = "PE")
	lines!(ax1, t_dns, vec(∫S✶zdz_dns), label = "BPE", linestyle = :dash)
	axislegend(ax1, position = :rb)
	ax2 = Axis(fig2[2, 1], title = "Available potential energy")
	lines!(ax2, t_dns, vec(∫Szdz_dns) .- vec(∫S✶zdz_dns), label = "APE")
	axislegend(ax2, position = :rt)
	fig2
end

# ╔═╡ d552957f-3c0d-40f2-b356-1914449c0e2d
let
	fig2 = Figure(size = (500, 1000))
	ax1 = Axis(fig2[1, 1], title = "Time change salinity potential energies")
	lines!(ax1, t_dns[2:end], dₜ∫Szdz_dns, label = "dₜPE")
	lines!(ax1, t_dns[2:end], dₜ∫S✶zdz_dns, label = "dₜBPE", linestyle = :dash)
	axislegend(ax1, position = :rb)
	ax2 = Axis(fig2[2, 1], title = "Time change available potential energy")
	lines!(ax2, t_dns[2:end], dₜ∫Szdz_dns .- dₜ∫S✶zdz_dns, label = "dₜAPE", color = :red)
	axislegend(ax2, position = :rt)
	fig2
end

# ╔═╡ 998bbe9a-5846-47a5-8c52-ab550f0c32fd
md"""
### Flux and BPE
"""

# ╔═╡ 1b4b28cb-5220-4a4f-ab98-7d19de52343c
let
	fig = Figure(size = (500, 1000))
	ax = Axis(fig[1, 1], title = "Salt flux and BPE")
	lines!(ax, t_dns[2:end], dₜ∫S✶zdz_dns, label = "dₜ∫S✶zdz")
	lines!(ax, t_dns[2:end], ∫dₜ∫Sdzdz_dns, label = "∫dₜ∫Sdzdz", linestyle = :dash)
	axislegend(ax, position = :rb)
	ax2 = Axis(fig[2, 1],
				title = "Absolute error between salinity flux and background salt energy",
				xlabel = "time (s)",
				ylabel = "Absolute error")
	lines!(ax2, t_dns[2:end], abs.(∫dₜ∫Sdzdz_dns - dₜ∫S✶zdz_dns))
	fig
end

# ╔═╡ cfb85301-23bc-4074-9a28-7b0a4952142f
md"""
## Density
"""

# ╔═╡ 218f3e2f-a448-4c8c-9101-e542ae8478c2
begin
	# Potential energy no gravity
	∫σzdz_dns = 9.81 * sum(σ₀_dns .* z_dns * Δz_dns, dims = 1)
	dₜ∫σzdz_dns = vec(diff(∫σzdz_dns, dims = 2)) ./ Δt_dns
 
	# Sorted density, density contet, density flux
	σ✶_dns = similar(σ₀_dns)
	∫σdz_dns = similar(σ₀_dns)
	for i ∈ eachindex(t_dns)
		σ✶_dns[:, i] = sort(σ₀_dns[:, i], rev = true)
		∫σdz_dns[:, i] = cumsum(reverse(σ✶_dns[:, i]) * Δz)
	end

	dₜ∫σdz_dns = diff(∫σdz_dns, dims = 2)
	∫dₜ∫σdzdz_dns = (vec(sum(dₜ∫σdz_dns * Δz_dns, dims = 1)) ./ Δt_dns) ./ 100

	# Background potential energy
	∫σ✶zdz_dns = 9.81 * sum(σ✶_dns .* z_dns * Δz_dns, dims = 1)
	dₜ∫σ✶zdz_dns = vec(diff(∫σ✶zdz_dns, dims = 2)) ./ Δt_dns

	nothing
end

# ╔═╡ 9833d0cc-d423-4901-abb3-bad2ad8fe895
md"""
### Time evolution
"""

# ╔═╡ 7718a697-809a-43af-9a9a-87da66cc655d
let
	fig, ax, hm = heatmap(t_dns, z_dns_plot, σ₀_dns', colormap = :dense)
	ax.title = "Density hovmoller"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm, label = "σ₀ (kg/m^3 )")
	fig
end

# ╔═╡ c0ceb3b3-0c79-4a69-ade1-8c6281b3bab1
begin
	t_slider_dns_2 = @bind timestep_dns_2 PlutoUI.Slider(eachindex(t_dns))
	md"""
	Slider for plotting density profiles (profile and sorted profile)
	$(t_slider_dns_2)
	"""
end

# ╔═╡ 68b20f1d-649b-4cb7-9d80-7f66e08691f4
let
	fig, ax = lines(σ₀_dns[:, timestep_dns_2], z_dns_plot, label = "Initial profile")
	lines!(ax, σ✶_dns[:, timestep_dns_2], z_dns_plot, label = "Sorted profile", linestyle = :dash)
	ax.title = "Density profiles t = $(t[timestep_dns_2] / 60) minutes"
	ax.xlabel = "Density (σ₀, kg/m^3)"
	ax.ylabel = "z (m)"
	axislegend(ax, position = :lb)
	fig
end

# ╔═╡ 551f307c-f4bc-4b3b-9704-e721e2429247
md"""
### Energetics
"""

# ╔═╡ 17615a2c-41d5-454b-b786-1c87191ffc31
let
	fig2 = Figure(size = (500, 1000))
	ax1 = Axis(fig2[1, 1], title = "Potential energies")
	lines!(ax1, t_dns, vec(∫σzdz_dns), label = "PE")
	lines!(ax1, t_dns, vec(∫σ✶zdz_dns), label = "BPE", linestyle = :dash)
	axislegend(ax1, position = :rb)
	ax2 = Axis(fig2[2, 1], title = "Available potential energy")
	lines!(ax2, t_dns, vec(∫σzdz_dns) .- vec(∫σ✶zdz_dns), label = "APE")
	axislegend(ax2, position = :rt)
	fig2
end

# ╔═╡ cb5cc791-f7e1-450e-bbc9-e0edb2b22498
let
	fig2 = Figure(size = (500, 1000))
	ax1 = Axis(fig2[1, 1], title = "Time change potential energies")
	lines!(ax1, t_dns[2:end], dₜ∫σzdz_dns, label = "dₜPE")
	lines!(ax1, t_dns[2:end], dₜ∫σ✶zdz_dns, label = "dₜBPE", linestyle = :dash)
	axislegend(ax1, position = :rb)
	ax2 = Axis(fig2[2, 1], title = "Time change available potential energy")
	lines!(ax2, t_dns[2:end], dₜ∫σzdz_dns .- dₜ∫σ✶zdz_dns, label = "dₜAPE", color = :red)
	axislegend(ax2, position = :rt)
	fig2
end

# ╔═╡ cb227724-0437-444d-98bc-be2f658c2702
md"""
### Flux and BPE
"""

# ╔═╡ 5e993f26-9ce3-4f11-83a8-6c6f2a78fc12
let
	fig = Figure(size = (500, 1000))
	ax = Axis(fig[1, 1], title = "Density flux and BPE")
	lines!(ax, t_dns[2:end], dₜ∫σ✶zdz_dns, label = "dₜ∫σ✶zdz")
	lines!(ax, t_dns[2:end], ∫dₜ∫σdzdz_dns, label = "∫dₜ∫Sdzdz", linestyle = :dash)
	axislegend(ax, position = :rb)
	ax2 = Axis(fig[2, 1],
				title = "Absolute error between density flux and background potential energy",
				xlabel = "time (s)",
				ylabel = "Absolute error (log10)")
	lines!(ax2, t_dns[2:end], log10.(abs.(∫dₜ∫σdzdz_dns - dₜ∫σ✶zdz_dns)))
	fig
end

# ╔═╡ cd5de5b1-33f9-4a32-8fc6-6724ce1570a8
md"""
## Diffusivity

### Salt or (temperature)
"""

# ╔═╡ f3b225dc-8cc6-444d-b6ae-7ef7798f3303
begin
	dSdz_dns = diff(reverse(S✶_dns, dims = 1), dims = 1) ./ Δz_dns
	replace!(dSdz_dns, 0 => NaN)
	κₛ_dns = dₜ∫Sdz_dns[2:end, :] ./ dSdz_dns[:, 2:end]
	zrange_dns = 690:710
	∫κₛ_dns = sum(κₛ_dns[.!isnan.(κₛ_dns)] * Δz_dns, dims = 1)
	κₛ_dns_mean = mean(κₛ_dns[.!isnan.(κₛ_dns)])
	nothing
end

# ╔═╡ 65c16680-0150-44c8-8cad-5a0023abf916
begin
	z_slider_dns = @bind z_dns_ PlutoUI.Slider(eachindex(z_dns[1:end-1]))
	md"""
	Choose a depth level to display the diffusivity calculated $(z_slider_dns).
	"""
end

# ╔═╡ 2bee25dc-8e7d-4476-9b13-b3b3760533bd
let
	# zrange = 680:720
	# fig, ax, hm = heatmap(t[2:end], z_dns[:], κₛ_dns[:, :]')
	# Colorbar(fig[1, 2], hm)
	# fig
	# fig, ax = series(t[20:end-1], κₛ[zrange, 20:end], labels = ["z = -$(round(i, digits = 2))m" for i ∈ z[zrange]], solid_color=:black)
	# ax.title = "Diffusivity estimates for salinity about the model interface"
	# ax.xlabel = "time (s)"
	# ax.ylabel = "κ (m2/s)"
	# axislegend(ax)
	# Legend(fig[1, 2], ax)
	fig, ax = lines(t_dns[20:end-1], κₛ_dns[z_dns_, 20:end], label = "z = -$(round(z_dns[z_dns_], digits = 2))", solid_color=:black)
	ax.title = "Diffusivity estimates for salinity about the model interface"
	ax.xlabel = "time (s)"
	ax.ylabel = "κ (m2/s)"
	hlines!(1e-7, label = "κ molecular", color = :red, linestyle = :dash)
	axislegend(ax)
	fig
end

# ╔═╡ 666c8467-6460-4ae9-adca-27c241ef3fdd
TableOfContents()

# ╔═╡ Cell order:
# ╟─a290de4c-f791-11ee-0e04-11e4e34b7428
# ╟─ac639feb-9ee4-43f4-acd9-466fe3478d40
# ╟─edd033e9-55ca-4b68-bf03-330baaa35e67
# ╟─d2a1ec90-16ef-47b2-9adb-8829681ad5d9
# ╟─0efe7042-82d9-44aa-be55-eb3634425bd5
# ╟─daf8bad7-444b-4da7-968b-32f29f1b7eba
# ╟─06bd6b0a-8d01-43c1-8560-f3c20972fe9e
# ╟─08a90e73-8312-4cd2-9c5d-8ff9829962bc
# ╟─ebff4d86-8013-4473-ae43-0fc2b4fe2015
# ╟─1b7fba62-557e-4c53-8f67-b45b16548e1f
# ╟─944526a6-5154-44ac-af06-c7dba206007c
# ╟─3e3eb9fd-8208-4f99-8cbb-730e4a501772
# ╟─ab6e5cd8-9a90-4cc2-8e83-330aa7ae8b30
# ╟─16562398-cd61-4697-b94a-931a5d44a2f3
# ╟─b3fc213d-03c1-4fb5-9c13-f3f6e5f1f648
# ╟─385c06a8-f60b-4f5b-8bd2-e2a56447397c
# ╟─2a164d51-759d-4341-811d-d1fd406c0c3d
# ╟─6cce82f4-b6a8-4be4-8e15-e93908d72eee
# ╟─febc3ebe-8d07-4c7d-aaec-fd04bb7378a2
# ╟─c98ca6ac-0af9-4bbf-88d7-8ef72a58b3e0
# ╟─4fea18d5-9fb5-40ff-8ea0-7e9bbf3cb9e1
# ╟─3217dc82-a988-448f-9bd4-ef5f62b75630
# ╟─e01432d2-627c-44cf-aba7-67d2094092d3
# ╟─1a9eee8b-545f-4da3-a931-49d447521e3c
# ╟─0121a899-a1bf-4fc8-9aa2-ce5c801753ec
# ╟─9a82d299-0274-4e68-9c4b-da1350e52fe1
# ╟─2642bc19-9c7c-45f5-b199-fe318da98101
# ╟─52a15ea8-600d-4bbf-8dee-def4895a4ded
# ╟─cd11a001-c306-4e16-8df5-8a51eb0a46ee
# ╟─4f5f7717-53fc-4252-804a-c0a3293df8d5
# ╟─0602c590-4435-4963-acd3-15798c9de894
# ╟─ebf2aa23-5cdd-4de7-a285-ce1534666899
# ╟─7a2b5232-e081-45ae-95be-7d219aae91ed
# ╟─ffdd23ec-190f-4307-adc9-7058bfcc2593
# ╟─b4dcdf1e-09ce-41e5-83de-90356dce3eeb
# ╟─7579eacb-fddd-40bb-a50f-a8f98cbbdbcd
# ╟─cf9aba94-ef48-4478-a49a-92e9da84e9d9
# ╟─d552957f-3c0d-40f2-b356-1914449c0e2d
# ╟─998bbe9a-5846-47a5-8c52-ab550f0c32fd
# ╟─1b4b28cb-5220-4a4f-ab98-7d19de52343c
# ╟─cfb85301-23bc-4074-9a28-7b0a4952142f
# ╟─218f3e2f-a448-4c8c-9101-e542ae8478c2
# ╟─9833d0cc-d423-4901-abb3-bad2ad8fe895
# ╟─7718a697-809a-43af-9a9a-87da66cc655d
# ╟─c0ceb3b3-0c79-4a69-ade1-8c6281b3bab1
# ╟─68b20f1d-649b-4cb7-9d80-7f66e08691f4
# ╟─551f307c-f4bc-4b3b-9704-e721e2429247
# ╟─17615a2c-41d5-454b-b786-1c87191ffc31
# ╟─cb5cc791-f7e1-450e-bbc9-e0edb2b22498
# ╟─cb227724-0437-444d-98bc-be2f658c2702
# ╟─5e993f26-9ce3-4f11-83a8-6c6f2a78fc12
# ╟─cd5de5b1-33f9-4a32-8fc6-6724ce1570a8
# ╟─f3b225dc-8cc6-444d-b6ae-7ef7798f3303
# ╟─65c16680-0150-44c8-8cad-5a0023abf916
# ╟─2bee25dc-8e7d-4476-9b13-b3b3760533bd
# ╟─666c8467-6460-4ae9-adca-27c241ef3fdd
