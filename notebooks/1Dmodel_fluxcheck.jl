### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ a290de4c-f791-11ee-0e04-11e4e34b7428
begin
	using Pkg
	Pkg.activate("..")
	using JLD2, Oceananigans, CairoMakie
end

# ╔═╡ ac639feb-9ee4-43f4-acd9-466fe3478d40
md"""
# Flux and background state check with 1D model

While computing the background state and the buoyancy flux for the model there were some discrepancies.
In this notebook I run equivalent 1D model and analysis to try and find out what might be going on with some of these calculations.

The main two experiments I will compare are the isothermal experiment and the cabbeling experiment (this is all I have done for DNS).
"""

# ╔═╡ edd033e9-55ca-4b68-bf03-330baaa35e67
md"""
## 1D model - isothermal

This isothermal run of the 1D model has the same setup as the full DNS experiment.
"""

# ╔═╡ d2a1ec90-16ef-47b2-9adb-8829681ad5d9
begin
	iso_data = "../1DModel/OneDModelOutput_cabbeling_cd10.jld2"
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
	dₜ∫Szdz = vec(diff(∫Szdz, dims = 2)) #./ Δt

	# Salt sorted, content and flux
	S_sorted = similar(S.data[1, 1, :, :])
	S_content = similar(S.data[1, 1, :, :])
	for i ∈ eachindex(t)
		S_sorted[:, i] = sort(interior(S, 1, 1, :, i), rev = true)	
		S_content[:, i] = cumsum(S_sorted[:, i] .* -Δz) # negative to sign match above
	end

	S_flux = diff(S_content, dims = 2) #./ Δt
	S_flux_int = vec(sum(S_flux * Δz, dims = 1))
	
	# Background potential energy
	∫S✶zdz = sum(S_sorted .* z * Δz, dims = 1)
	dₜ∫S✶zdz = vec(diff(∫S✶zdz, dims = 2)) #./ Δt
	
	nothing
end

# ╔═╡ cc7d7952-eb79-47df-bd0a-9ba22963766a
S_sorted

# ╔═╡ 38d0188e-d6d0-4d44-bb75-638f28f63aaf
interior(S, 1, 1, :, :)

# ╔═╡ 0efe7042-82d9-44aa-be55-eb3634425bd5
let
	fig, ax, hm = heatmap(t, z, interior(σ₀, 1, 1, :, :)', colormap = :dense)
	ax.title = "Density hovmoller"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm, label = "σ₀ (kg/m3 )")
	fig
end

# ╔═╡ ebff4d86-8013-4473-ae43-0fc2b4fe2015
md"""
## Energetics

Plot of potential, background potential and available potential energies using only salt tracer.
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
	lines!(ax2, t[2:end], dₜ∫S✶zdz, label = "BPE")
	lines!(ax2, t[2:end], dₜ∫Szdz .- dₜ∫S✶zdz, label = "APE")
	ax2.title = "Potential energies"
	axislegend(ax2, position = :rb)
	fig2
end

# ╔═╡ 3e3eb9fd-8208-4f99-8cbb-730e4a501772
md"""
## Background state and diffusive salt flux

Changes to the background potential energy can **only occur due to vertical diffusive flux** so we should have equality between
```math
\frac{\mathrm{d}}{\mathrm{d}t}E_{b} \quad \mathrm{and} \quad \frac{\mathrm{d}}{\mathrm{d}t}\iint S^{*} dz dz
```
"""

# ╔═╡ ab6e5cd8-9a90-4cc2-8e83-330aa7ae8b30
let
	fig, ax = lines(t[2:end], dₜ∫S✶zdz, label = "BPE")
	lines!(ax, t[2:end], S_flux_int, label = "Salt flux")
	ax.title = "Salt flux and BPE"
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ e14dbf9d-3a92-4051-b99f-8a7c9eb7de5f
abs.(S_flux_int - dₜ∫S✶zdz)

# ╔═╡ 8d029027-9e35-40cc-a7c7-dbbb94878724
md"""
## Other things
"""

# ╔═╡ df8c87bc-ba1d-402d-8d40-089ba0f7447f
let
	fig, ax, hm = heatmap(Matrix(S_flux[:, :])')
	Colorbar(fig[1, 2], hm)
	fig
end

# ╔═╡ 42f95fcf-be3c-4c8d-8226-5dd33bf73b7a
begin
	S_diff = S_flux[2:end, :] ./ diff(S_sorted[:, 2:end], dims = 1)
	replace!(S_diff, Inf => 0)
	replace!(S_diff, -Inf => 0)
end

# ╔═╡ c5b279df-1031-433f-acf4-9c32c07fcfcb
let
	# fig, ax, hm = heatmap(Matrix(S_diff[695:715, :])')
	# Colorbar(fig[1, 2], hm)
	# fig
	lines(vec(S_diff[49, :]))
end

# ╔═╡ Cell order:
# ╟─a290de4c-f791-11ee-0e04-11e4e34b7428
# ╟─ac639feb-9ee4-43f4-acd9-466fe3478d40
# ╟─edd033e9-55ca-4b68-bf03-330baaa35e67
# ╠═d2a1ec90-16ef-47b2-9adb-8829681ad5d9
# ╠═cc7d7952-eb79-47df-bd0a-9ba22963766a
# ╠═38d0188e-d6d0-4d44-bb75-638f28f63aaf
# ╟─0efe7042-82d9-44aa-be55-eb3634425bd5
# ╟─ebff4d86-8013-4473-ae43-0fc2b4fe2015
# ╟─944526a6-5154-44ac-af06-c7dba206007c
# ╟─3e3eb9fd-8208-4f99-8cbb-730e4a501772
# ╠═ab6e5cd8-9a90-4cc2-8e83-330aa7ae8b30
# ╠═e14dbf9d-3a92-4051-b99f-8a7c9eb7de5f
# ╟─8d029027-9e35-40cc-a7c7-dbbb94878724
# ╟─df8c87bc-ba1d-402d-8d40-089ba0f7447f
# ╠═42f95fcf-be3c-4c8d-8226-5dd33bf73b7a
# ╠═c5b279df-1031-433f-acf4-9c32c07fcfcb
