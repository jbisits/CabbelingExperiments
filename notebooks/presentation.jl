### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 46063034-7f6c-11ee-18e4-7bb553b5bb66
begin 
	using Pkg
	Pkg.activate("..")
	using CairoMakie, PlutoUI, GibbsSeaWater
end

# ╔═╡ 8ce9ed4d-0a8f-410c-a021-91125d65a9a5
md"""
# The cabbeling instability in Direct Numerical Simulations
"""

# ╔═╡ 8d54f303-4b14-490d-bdb3-aa88c3ccf7bd
md"""
# Non linearities in the equation of state - 1

- Non-linearities in equation of state impact the formation of sea-ice and the internal pycnocline in the Southern Ocean. **References**.
- Whilst not leading order, cabbeling and thermobaricity are two important non-linear processes influencing the formation of intermediate waters and the layering of the deep ocean respectively.
- Broad sentence why we keep studying them.

"""

# ╔═╡ b1f2293e-9a64-4f27-a07b-f2b1f6ab6611
begin
	S_star, Θ_star = 34.7, 0.5
	S₀ᵘ = 34.58
	Θᵘ = -1.5
	slope = (Θᵘ - Θ_star) / (S₀ᵘ - S_star)
	S_mix = range(S₀ᵘ, S_star, step = 0.000001)
	Θ_mix = @. Θᵘ + (slope) * (S_mix - S₀ᵘ)
	ρ_mix = gsw_rho.(S_mix, Θ_mix, 0)
	max_rho, max_rho_idx = findmax(ρ_mix)
	S_max, Θ_max = S_mix[max_rho_idx], Θ_mix[max_rho_idx]
	Δρ_mix = max_rho - gsw_rho(S_star, Θ_star, 0)

	N = 2000
	S_range, Θ_range = range(34.52, 34.72, length = N), range(-2, 1, length = N)
	S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
	ρ = gsw_rho.(S_grid, Θ_grid, 0)
	ρ_star = gsw_rho(S_star, Θ_star, 0)
	ρ_s = gsw_rho(S₀ᵘ, Θᵘ, 0)
	find_Θ = findfirst(Θ_range .> -1.5)
	find_S = findfirst(ρ[find_Θ, :] .> ρ_star)
	S_iso, Θ_iso = S_range[find_S], Θ_range[find_Θ]
	gsw_rho(S_iso, Θ_iso, 0)
	αₗ, βₗ = gsw_alpha(S_star, Θ_star, 0), gsw_beta(S_star, Θ_star, 0)
	m_initial = βₗ / αₗ
	Θ_linear_initial = @. Θ_star + m_initial * (S_range - S_star)
	αₘ, βₘ = gsw_alpha(S_max, Θ_max, 0), gsw_beta(S_max, Θ_max, 0)
	m = βₘ / αₘ
	Θ_linear = @. Θ_max + m * (S_range - S_max)
	TSfig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(TSfig[1, 1];
			  title = "Cabbeling effect in S-Θ space",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], color = :red, linewidth = 0.8, labelsize = 18, linestyle = :dot)
	scatter!(ax, [S_star], [Θ_star], color = :red, label = "Deep water")
	scatter!(ax, [S₀ᵘ], [Θᵘ], color = :blue, label = "Shallow water")
	lines!(ax, S_mix, Θ_mix, color = :purple, linestyle = :dot, label = "Mixed water")
	axislegend(ax, position = :lt)
	md"""
	# Cabbeling - 2
	
	- The dominant non-linear process in the upper ocean (to around ``z = 1000~m``) is *cabbeling*.
	- Cabbeling is the gain in density that occurs when water masses of differing temperature and salinity but equal density mix.

	$(TSfig)
	"""
end

# ╔═╡ 57bf32e9-c9fc-444b-8ef9-d5fdd8ecf4e4
md"""
# Cabbeling - 3

- In my first project we developed a method to diagnose when a profile will be unstable to cabbeling,
```math
\Delta\rho > \rho\left(S^{*} + \frac{\alpha^{*}}{\beta^{*}}\Delta\Theta, \Theta^{*} + \Delta\Theta, \overline{p} \right) - \rho\left(S^{*}, \Theta^{*}, \overline{p}\right).
```
$(LocalResource("../../../Papers/PhD-paper1-CabbelingInstability/Plots_v4/fig6_ECCOpdfs.png"))
"""

# ╔═╡ 5987f662-879f-4394-aac9-a45634bf8ab3
md"""
# Cabbeling - 4

- This work provided evidence that cabbeling is influencing *how unstable temperature inverted profiles can become prior to an instability forming*... but
- we still could not say anything about the mixing that was taking place.
"""

# ╔═╡ 5bc2dbc2-f295-4389-be7a-b172cfbe6702
md"""
# Investigating the cabbeling instability in Direct Numerical Simulations

- Direct Numerical Simulations (DNS) simulate all scales on which turbulence can occur.
- By keeping key non-dimensional parameters close to realistic values we can simulate the mixing that is taking place at very small scales.
- To our knowledge no one has used DNS to explore cabbeling.
"""

# ╔═╡ 1b10d122-b550-47b8-aa78-1c073d14167c
md"""
# Simulating turbulence on all scales
- Eddies form on many different length scales with larger eddies creating smaller eddies and transferring energy to them.
- The transfer of energy from largest to smallest scale is called the energy cascade.
- At all length scales, energy is being **dissapated by viscosity.**
- It follows that there is a length scale at which viscous effects dominate and eddies will not form,
```math
\eta = \left(\frac{\nu^{3}}{\epsilon}\right)^{1/4} \tag{1}
```
- ``(1)`` is known as the Kolmogorov length scale and represents the smallest length of an eddy in a turbulent flow.
- In DNS the grid resolution must satisfy ``(1)`` at all points in space and time to be simulating turbulence on all length scales. 
"""

# ╔═╡ dfd413ab-4184-4d4b-a044-df0802d8647d
md"""
# Non-dimensional parameters
- DNS requires extremely high resolution to be at the Kolmogorov scale meaning other aspects of the flow, namely viscosity, may be altered to limit computational demand.
- Increasing the magnitude of the viscosity and adjusting other model values so that non-dimensional parameters remain the same, ensures the characteristics of the flow are not altered.
"""

# ╔═╡ 39ba4eb8-4040-480a-9ea5-08d538f1f2ea
md"""
# Non-dimensional parameters
- DNS requires extremely high resolution to be at the Kolmogorov scale meaning other aspects of the flow, namely viscosity, may be altered to limit computational demand.
- Increasing the magnitude of the viscosity and adjusting other model values so that non-dimensional parameters remain the same, ensures the characteristics of the flow are not altered.
- Prandtl number: ``\mathrm{Pr} = \frac{\nu}{\kappa_{\Theta}}``
- Schmidt number: ``\mathrm{Sc} = \frac{\nu}{\kappa_{S}}``
- Diffusivity ratio: ``\kappa_{R} = \frac{\kappa_{\Theta}}{\kappa_{S}}``
"""

# ╔═╡ 4c7769d2-852f-4c66-ad6c-31f453aedda9
md"""
# What have we actually done!
"""

# ╔═╡ 9feac837-c922-41b3-9e14-ba704aa2386b
TableOfContents(title = "Slides")

# ╔═╡ Cell order:
# ╟─46063034-7f6c-11ee-18e4-7bb553b5bb66
# ╟─8ce9ed4d-0a8f-410c-a021-91125d65a9a5
# ╟─8d54f303-4b14-490d-bdb3-aa88c3ccf7bd
# ╟─b1f2293e-9a64-4f27-a07b-f2b1f6ab6611
# ╟─57bf32e9-c9fc-444b-8ef9-d5fdd8ecf4e4
# ╟─5987f662-879f-4394-aac9-a45634bf8ab3
# ╟─5bc2dbc2-f295-4389-be7a-b172cfbe6702
# ╟─1b10d122-b550-47b8-aa78-1c073d14167c
# ╟─dfd413ab-4184-4d4b-a044-df0802d8647d
# ╟─39ba4eb8-4040-480a-9ea5-08d538f1f2ea
# ╠═4c7769d2-852f-4c66-ad6c-31f453aedda9
# ╟─9feac837-c922-41b3-9e14-ba704aa2386b
