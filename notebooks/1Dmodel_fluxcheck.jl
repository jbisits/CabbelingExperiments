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
	using JLD2, Oceananigans, CairoMakie, PlutoUI
end

# ╔═╡ ac639feb-9ee4-43f4-acd9-466fe3478d40
md"""
# Flux and background state check with 1D model

While computing the background state and the buoyancy flux for the model there were some discrepancies.
In this notebook I run equivalent 1D model and analysis to try and find out what might be going on with some of these calculations.

The main two experiments I will compare are the isothermal experiment and the cabbeling experiment (this is all I have done for DNS).

## Model runs

I have three model runs I am comparing:

1. isothermal, uniform ``T = 0.5^{\circ}C`` (there is also a lower resolution one of these experiments)
2. cabbeling with convective diffusivity ``\kappa_{c} = 1\mathrm{m}^{2}\mathrm{s}^{-1}``
2. cabbeling with convective diffusivity ``\kappa_{c} = 10\mathrm{m}^{2}\mathrm{s}^{-1}``

In the isothermal case, the background state should be equal to the profile so the BPE and the PE should be equal (or very close to equal) throughout the length of the because the mixed water at the interface will not be denser/lighter than the lower/upper layers.

In the cabbeling case, the background state and profile should start roughly equal but after the mixing has created some denser water the PE and the BPE should separate until when/if lower layer has all been transformed to the maximum density.
"""

# ╔═╡ edd033e9-55ca-4b68-bf03-330baaa35e67
md"""
## Salinity
"""

# ╔═╡ 0f390348-56c2-46ee-98c3-0cd7fc09b1d3
@bind experiment Select(["isothermal", "cabbeling_cd1", "cabbeling_cd10"])

# ╔═╡ d2a1ec90-16ef-47b2-9adb-8829681ad5d9
begin
	iso_data = "../1DModel/OneDModelOutput_"*experiment*".jld2"
	S = FieldTimeSeries(iso_data, "S")
	T = FieldTimeSeries(iso_data, "T")
	σ₀ = FieldTimeSeries(iso_data, "σ")
	t = S.times
	Δt = diff(t)
	z = znodes(S)
	Δz = S.grid.Δzᵃᵃᶜ
	
	#Computations to check

	# Potential energy for salt
	∫Szdz = sum(interior(S, 1, 1, :, :) .* z * Δz, dims = 1)
	dₜ∫Szdz = vec(diff(∫Szdz, dims = 2)) ./ Δt

	# Salt sorted, content and flux
	S_sorted = similar(S.data[1, 1, :, :])
	S_content = similar(S.data[1, 1, :, :])
	for i ∈ eachindex(t)
		S_sorted[:, i] = sort(interior(S, 1, 1, :, i), rev = true)	
		S_content[:, i] = -cumsum(S_sorted[:, i] * Δz) # negative to sign match above
	end

	# S_flux = diff(S_content, dims = 2)
	# S_flux_int = vec(sum(S_flux * Δz, dims = 1)) #./ Δt
	S_flux = vec(sum(S_content * Δz, dims = 1))
	S_flux_int = diff(S_flux) ./ Δt
	
	# Background potential energy
	∫S✶zdz = sum(S_sorted .* z * Δz, dims = 1)
	dₜ∫S✶zdz = vec(diff(∫S✶zdz, dims = 2)) ./ Δt
	
	nothing
end

# ╔═╡ 0efe7042-82d9-44aa-be55-eb3634425bd5
let
	fig, ax, hm = heatmap(t, z, interior(S, 1, 1, :, :)', colormap = :haline)
	ax.title = "Salinity hovmoller"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm, label = "S (g/kg )")
	fig
end

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
	fig, ax = lines(interior(S, 1, 1, :, timestep), z, label = "Initial profile")
	lines!(ax, S_sorted[:, timestep], z, label = "Sorted profile", linestyle = :dash)
	ax.title = "Salinity profiles t = $(t[timestep] / 60) minutes"
	ax.xlabel = "Salinity (g/kg)"
	ax.ylabel = "z (m)"
	axislegend(ax)
	fig
end

# ╔═╡ ebff4d86-8013-4473-ae43-0fc2b4fe2015
md"""
### Energetics using salinity only

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

# ╔═╡ 944526a6-5154-44ac-af06-c7dba206007c
begin
	fig2, ax2 = lines(t[2:end], dₜ∫Szdz, label = "PE")
	lines!(ax2, t[2:end], dₜ∫S✶zdz, label = "BPE", linestyle = :dash)
	lines!(ax2, t[2:end], dₜ∫Szdz .- dₜ∫S✶zdz, label = "APE", linestyle = :dot)
	ax2.title = "Potential energies"
	axislegend(ax2, position = :rb)
	fig2
end

# ╔═╡ 3e3eb9fd-8208-4f99-8cbb-730e4a501772
md"""
### Background salinity state and diffusive salt flux

Changes to the background potential energy can **only occur due to vertical diffusive flux** so we should have equality between
```math
\frac{\mathrm{d}}{\mathrm{d}t}E_{b} \quad \mathrm{and} \quad \frac{\mathrm{d}}{\mathrm{d}t}\iint S^{*} dz dz
```
"""

# ╔═╡ ab6e5cd8-9a90-4cc2-8e83-330aa7ae8b30
let
	fig, ax = lines(t[2:end], dₜ∫S✶zdz, label = "Salt BPE")
	lines!(ax, t[2:end], S_flux_int, label = "Salt flux", linestyle = :dash)
	ax.title = "Salt flux and BPE"
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ e14dbf9d-3a92-4051-b99f-8a7c9eb7de5f
let
	fig, ax = lines(t[2:end], abs.(S_flux_int - dₜ∫S✶zdz))
	ax.xlabel = "time (s)"
	ax.ylabel = "Absolute error"
	ax.title = "Absolute error between salinity flux and background salt energy"
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
	g = 1 #9.81
	∫σzdz = g * sum(interior(σ₀, 1, 1, :, :) .* z * Δz, dims = 1)
	dₜ∫σzdz = vec(diff(∫σzdz, dims = 2)) ./ Δt

	# Sorted density, density contet, density flux
	σ_sorted = similar(σ₀.data[1, 1, :, :])
	σ_content = similar(σ₀.data[1, 1, :, :])
	for i ∈ eachindex(t)
		σ_sorted[:, i] = sort(interior(σ₀, 1, 1, :, i), rev = true)	
		σ_content[:, i] = cumsum(σ_sorted[:, i] * -Δz) # negative to sign match above
	end

	σ_flux = diff(σ_content, dims = 2) 
	σ_flux_int = vec(sum(σ_flux * Δz, dims = 1)) ./ Δt
	
	# Background potential energy
	∫σ✶zdz = g * sum(σ_sorted .* z * Δz, dims = 1)
	dₜ∫σ✶zdz = vec(diff(∫σ✶zdz, dims = 2)) ./ Δt

	nothing
end

# ╔═╡ 2a164d51-759d-4341-811d-d1fd406c0c3d
let
	fig, ax, hm = heatmap(t, z, interior(σ₀, 1, 1, :, :)', colormap = :dense)
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
	fig, ax = lines(interior(σ₀, 1, 1, :, timestep2), z, label = "Initial profile")
	lines!(ax, σ_sorted[:, timestep2], z, label = "Sorted profile", linestyle = :dash)
	ax.title = "Density profiles t = $(t[timestep2] / 60) minutes"
	ax.xlabel = "Density (σ₀, kg/m^3)"
	ax.ylabel = "z (m)"
	axislegend(ax)
	fig
end

# ╔═╡ 3217dc82-a988-448f-9bd4-ef5f62b75630
let
	fig2, ax2 = lines(t[2:end], dₜ∫σzdz, label = "PE")
	lines!(ax2, t[2:end], dₜ∫σ✶zdz, label = "BPE", linestyle = :dash)
	#lines!(ax2, t[2:end], dₜ∫σzdz .- dₜ∫σ✶zdz, label = "APE", linestyle = :dot)
	ax2.title = "Potential energies"
	axislegend(ax2, position = :rb)
	ax2.xlabel = "time (s)"
	ax2.ylabel = "pseudo-Watts"
	fig2
end

# ╔═╡ 1a9eee8b-545f-4da3-a931-49d447521e3c
let
	fig, ax = lines(t[2:end], dₜ∫σ✶zdz, label = "BPE")
	lines!(ax, t[2:end], σ_flux_int, label = "Density flux", linestyle = :dash)
	ax.title = "Density flux and BPE"
	ax.xlabel = "time (s)"
	ax.ylabel = "pseduo-Watts"
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ 725627af-aafc-42de-9462-67577e01de98
let
	fig, ax = lines(t[2:end], log10.(abs.(σ_flux_int - dₜ∫σ✶zdz)))
	ax.xlabel = "time (s)"
	ax.ylabel = "Absolute error (log10)"
	ax.title = "Absolute error between density flux and background energy"
	fig
end

# ╔═╡ 666c8467-6460-4ae9-adca-27c241ef3fdd
TableOfContents()

# ╔═╡ Cell order:
# ╟─a290de4c-f791-11ee-0e04-11e4e34b7428
# ╟─ac639feb-9ee4-43f4-acd9-466fe3478d40
# ╟─edd033e9-55ca-4b68-bf03-330baaa35e67
# ╟─0f390348-56c2-46ee-98c3-0cd7fc09b1d3
# ╟─d2a1ec90-16ef-47b2-9adb-8829681ad5d9
# ╟─0efe7042-82d9-44aa-be55-eb3634425bd5
# ╟─06bd6b0a-8d01-43c1-8560-f3c20972fe9e
# ╟─08a90e73-8312-4cd2-9c5d-8ff9829962bc
# ╟─ebff4d86-8013-4473-ae43-0fc2b4fe2015
# ╟─944526a6-5154-44ac-af06-c7dba206007c
# ╟─3e3eb9fd-8208-4f99-8cbb-730e4a501772
# ╟─ab6e5cd8-9a90-4cc2-8e83-330aa7ae8b30
# ╟─e14dbf9d-3a92-4051-b99f-8a7c9eb7de5f
# ╟─16562398-cd61-4697-b94a-931a5d44a2f3
# ╟─b3fc213d-03c1-4fb5-9c13-f3f6e5f1f648
# ╟─2a164d51-759d-4341-811d-d1fd406c0c3d
# ╟─6cce82f4-b6a8-4be4-8e15-e93908d72eee
# ╟─febc3ebe-8d07-4c7d-aaec-fd04bb7378a2
# ╟─3217dc82-a988-448f-9bd4-ef5f62b75630
# ╟─1a9eee8b-545f-4da3-a931-49d447521e3c
# ╟─725627af-aafc-42de-9462-67577e01de98
# ╟─666c8467-6460-4ae9-adca-27c241ef3fdd
