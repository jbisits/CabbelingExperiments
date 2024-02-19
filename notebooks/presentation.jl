### A Pluto.jl notebook ###
# v0.19.36

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
# Non linearities in the equation of state

- Non-linearities in equation of state impact the formation of sea-ice (Roquet et al. 2022) and the internal pycnocline (Klocker et al. 2023) in the Southern Ocean. 
- Whilst not leading order, cabbeling and thermobaricity are two important non-linear processes influencing the formation of intermediate waters and the layering of the deep ocean respectively (Nycander et al. 2015).
- How these non-linear processes influence the vertical exchanges of heat and carbon is important for understanding how the ocean helps regulate the global climate.
"""

# ╔═╡ b1f2293e-9a64-4f27-a07b-f2b1f6ab6611
let
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
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], color = :black, linewidth = 0.8, labelsize = 18, linestyle = :dot, label = "Isopycnal")
	scatter!(ax, [S_star], [Θ_star], color = :red, label = "Deep water")
	scatter!(ax, [S₀ᵘ], [Θᵘ], color = :blue, label = "Shallow water")
	lines!(ax, S_mix, Θ_mix, color = :purple, linestyle = :dash, label = "Mixed water")
	axislegend(ax, position = :lt)
	#save("STfig.png", TSfig)
	md"""
	# Cabbeling
	
	- The dominant non-linear process in the upper ocean (to around ``z = 1000~m``) is *cabbeling*.
	- Cabbeling is the gain in density that occurs when water parcels of differing temperature and salinity mix to form *denser water*.

	$(TSfig)
	"""
end

# ╔═╡ 57bf32e9-c9fc-444b-8ef9-d5fdd8ecf4e4
md"""
# Cabbeling

- In my first project we developed a method to diagnose when a profile will be unstable to cabbeling,
```math
\Delta\rho > \rho\left(S^{*} + \frac{\alpha^{*}}{\beta^{*}}\Delta\Theta, \Theta^{*} + \Delta\Theta, \overline{p} \right) - \rho\left(S^{*}, \Theta^{*}, \overline{p}\right).
```
$(LocalResource("../../../Papers/PhD-paper1-CabbelingInstability/fig6_ECCOpdfs.png"))
"""

# ╔═╡ 5987f662-879f-4394-aac9-a45634bf8ab3
let
	data = [vcat(fill(0, 5), fill(1, 5)) vcat(fill(0, 5), fill(1, 5))]
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1],
			  aspect = 1/3,
	          ylabel = "← z",
        	  yticksvisible = false,
        	  yticklabelsvisible = false,
        	  ygridvisible = false)
	hidexdecorations!(ax)
	hidespines!(ax)
	hm = heatmap!(ax, 1:2, 1:10, data'; colormap = [:red, :blue])
	text!(ax, 1.3, 7.5,  text = L"ρ^{-}", fontsize = 34)
	text!(ax, 1.3, 2.5,  text = L"ρ^{+}", fontsize = 34)
	ax2 = Axis(fig[1, 2],
			  aspect = 1/3,
	          ylabel = "← z",
        	  yticksvisible = false,
        	  yticklabelsvisible = false,
        	  ygridvisible = false)
	hidexdecorations!(ax2)
	hideydecorations!(ax2)
	hidespines!(ax2)
	data2 = [vcat(fill(0, 4), 1, fill(2, 4)) vcat(fill(0, 4), 1, fill(2, 4))]
	hm = heatmap!(ax2, 1:2, 1:10, data2'; colormap = [:red, :purple, :blue])
	text!(ax2, 1.3, 7.5,  text = L"ρ^{-}", fontsize = 34)
	text!(ax2, 1.3, 2.5,  text = L"ρ^{+}", fontsize = 34)
	text!(ax2, 1.1, 5.35,  text = L"ρ_{\mathrm{mixed}}", fontsize = 34)
	fig
	md"""
	# Cabbeling
	
	- This work provided evidence that cabbeling is influencing *how unstable temperature inverted profiles can become prior to an instability forming*... but
	- we still could not say anything about the mixing that was taking place at all scales and if cabbeling is driving any convective mixing.
	$(fig)
	"""
end

# ╔═╡ 5bc2dbc2-f295-4389-be7a-b172cfbe6702
md"""
# Investigating the cabbeling instability in Direct Numerical Simulations

- Direct Numerical Simulations (DNS) simulate all length scales on which turbulence can occur.
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

# ╔═╡ 39ba4eb8-4040-480a-9ea5-08d538f1f2ea
md"""
# Non-dimensional parameters
- DNS requires extremely high resolution to be at the Kolmogorov scale meaning other aspects of the flow, namely viscosity, may be altered to limit computational demand.
- Increasing the magnitude of the viscosity and adjusting other model values so that non-dimensional parameters remain the same, ensures the characteristics of the flow are not altered.
- Prandtl number: ``\mathrm{Pr} = \frac{\nu}{\kappa_{\Theta}}``
- Schmidt number: ``\mathrm{Sc} = \frac{\nu}{\kappa_{S}}``
- Diffusivity ratio: ``\kappa_{R} = \frac{\kappa_{\Theta}}{\kappa_{S}}``
- For ``\mathrm{Pr} > 1`` the Kolmogorov scale ``(1)`` is replaced by the Batchelor scale
```math
\lambda_{B} = \frac{\eta}{\sqrt{\mathrm{Pr}}}.
```
"""

# ╔═╡ ec83df92-8beb-4808-b1ac-8db61b81bae6
md"""
# Our simulation setup
- Using [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/dev/) we setup a Non-hydrostatic model that solves the incompressible Navier Stokes equations under the Boussinesq approximation.
- **Domain:** ``0.1m\times 0.1m\times 1m``.
- **Boundary conditions:** horizontally periodic, "zero-flux" vertically.
- **Kinematic viscosity:** ``1\times 10^{-6}m^{2}s^{-1}``.
- **Isotropic (molecular) diffusivity values:** ``\kappa_{S} = \kappa_{\Theta} = 1\times 10^{-7}m^{2}s^{-1}.``
- **Equation of State:** 55 term polynomial approximation to full EOS appropriate for Boussinesq models.
"""

# ╔═╡ e053dce3-1d3b-4012-afbe-65b1d754d8cb
md"""
# Model resolution
- To simulate turbulence on all scales occuring in our domain we need to have resolution less than the Kolmogorov scale.
- To achieve this we have **model resolution:** `Nx = 124, Ny = 124, Nz = 1400` to get uniform horizontal resolution ≈ 0.8mm and uniform vertical velocity ≈ 0.7mm.
- After running the simulations we found we had Kolmogorov length scale ≈ 1.8mm and Batchelor scale ≈ 0.6mm.
"""

# ╔═╡ 4c7769d2-852f-4c66-ad6c-31f453aedda9
let
	S_star, Θ_star = 34.7, 0.5
	S_isothermal = 34.69431424
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
			  title = "DNS initial conditions",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], color = :black, linewidth = 0.8, labelsize = 18, linestyle = :dot, label = "Isopycnal")
	scatter!(ax, [S_star], [Θ_star], color = :red, label = "Deep water")
	scatter!(ax, [S₀ᵘ], [Θᵘ], color = :blue, label = "Cabbeling shallow water")
	lines!(ax, S_mix, Θ_mix, color = :purple, linestyle = :dash, label = "Mixed water")
	scatter!(ax, [S_isothermal], [Θ_star], color = :yellow, label = "Isothermal shallow water")
	slope_iso = (Θ_star - Θ_star) / (S_isothermal - S_star)
	S_mix2 = range(S_isothermal, S_star, step = 0.000001)
	Θ_mix2 = @. Θ_star + (slope_iso) * (S_mix2 - S_isothermal)
	lines!(ax, S_mix2, Θ_mix2, color = :orange, linestyle = :dash, label = "Mixed water")
	axislegend(ax, position = :lt)
	md"""
	# Initial conditions
	- Two layer system that is statically stable.
	- Deeper layer constant across all experiments at ``S = 34.7`` gkg⁻¹, ``\Theta = 0.5`` °C
	- Small (order ``2\times 10^{-4}``) amount of random noise added about the interface of the two layer in salinity field to kick off mixing.
	- Test the case of two tracers vs one tracer mixing
	- **Cabbeling:** Cold/fresh shallow layer, ``S = 34.58`` gkg⁻¹, ``\Theta = -1.5`` °C
	- **Isothermal:** Fresh shallow layer, ``S_{\mathrm{shallow}} = 34.694`` gkg⁻¹
	$(TSfig)
	"""
	#save("DNS_ICS.png", TSfig)
end

# ╔═╡ e26daadf-0134-463a-a182-2a63e5c73276
md"""
# Isothermal experiment
$(LocalResource("../data/animations/Test cases/isothermal_density.mp4"))
"""

# ╔═╡ 0bffc2f9-1aa6-4cf6-9e3f-4010c97a7052
md"""
# Cabbeling experiment
$(LocalResource("../data/animations/Test cases/cabbeling_density.mp4"))
"""

# ╔═╡ 3dc7a488-12a3-41e7-a4bf-a715072c800d
begin
	shallow_salinity = @bind S₀ᵘ PlutoUI.Slider(range(34.551, 34.6, step = 0.001))
	nothing
end

# ╔═╡ ae8f6de0-7960-4beb-90af-7ea8e479cb56
let
	S_star, Θ_star = 34.7, 0.5
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
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1];
			  title = "Mixing water masses",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], color = :red, linewidth = 0.8, labelsize = 18, linestyle = :dot)
	lines!(ax, S_range, Θ_linear_initial, color = :red, linestyle = :dash, label = "Linear density at initial deep water", linewidth = 0.8)
	scatter!(ax, [S_star], [Θ_star], color = :red, label = "Deep water")
	scatter!(ax, [S₀ᵘ], [Θᵘ], color = :blue, label = "Shallow water")
	lines!(ax, S_mix, Θ_mix, color = :purple, linestyle = :dot, label = "Mixed water")
	scatter!(ax, S_mix[max_rho_idx], Θ_mix[max_rho_idx], color = :green, label = "Maximum ρ from mixing")
	contour!(ax, S_range, Θ_range, ρ', levels = [max_rho], color = :green, linestyle = :dot, linewidth = 0.8)
	lines!(ax, S_range, Θ_linear, color = :green, linewidth = 0.8, linestyle = :dash, label = "Linear density at new deep water")
	axislegend(ax, position = :lt)
	fig

	md"""
	# Maximum density from mixing
	- Set the shallow water at ``\Theta = `` $(Θᵘ)°C and vary ``S`` between stable to cabbeling and isohaline with the slider

	``S`` = $(shallow_salinity) = $(S₀ᵘ)gkg⁻¹.

	- The *maximum density after mixing* is $(round(max_rho, digits = 5)) kgm⁻³ which is a gain of $(round(Δρ_mix, digits = 5))kgm⁻³.
	- The new maximum density is at salinity $(round(S_mix[max_rho_idx], digits = 3))gkg⁻¹ and temperature $(round(Θ_mix[max_rho_idx], digits = 2))°C.

	- Slope of mixing line is $(slope) and linearised density slope is $(m).

	$(fig)
	"""


end

# ╔═╡ e46fcd1b-ccce-4089-8b49-66e93661fdc1
md"""
# New stable state

- Provided there is enough shallow water to bring all the deep water to the maximum density a new stable state will be reached where the cabbeling instability has been mixed away.
- The relative amounts of shallow and deep water needed to achieve a new stable to cabbeling state can be found by solving the linear system
```math
\begin{aligned}
a\begin{pmatrix}S_\mathrm{shallow} \\ \Theta_{\mathrm{shallow}} \end{pmatrix} + b\begin{pmatrix} S_\mathrm{deep} \\ \Theta_{\mathrm{deep}} \end{pmatrix} &= \begin{pmatrix} S_{\mathrm{mix}} \\ \Theta_{\mathrm{mix}} \end{pmatrix} \\
\begin{pmatrix}S_\mathrm{shallow} & S_{\mathrm{deep}} \\ \Theta_{\mathrm{shallow}} & \Theta_{\mathrm{deep}} \end{pmatrix} \begin{pmatrix} a \\ b \end{pmatrix} &= 
\begin{pmatrix} S_{\mathrm{mix}} \\ \Theta_{\mathrm{mix}} \end{pmatrix}.
\end{aligned}
```
- If there is not enough shallow water to bring all the deep water to the maximum density, the convective mixing will continue until there is no shallow water left.
"""

# ╔═╡ f82ebb94-c94d-42cc-bd16-92b123ba6498
md"""
# Longer experiments
- We test these hypotheses by running a longer experiment with an equal amount of shallow and deep water to see if the predicted maximum density is what the deeper water is transformed to and if the new stable to cabbeling state is reached.

$(LocalResource("../data/animations/cabbeling_600min/density.mp4"))
"""

# ╔═╡ bafe1643-5d0e-4d16-87ca-5edffef41c8c
md"""
# Tracers

$(LocalResource("../data/animations/cabbeling_600min/tracers_600min.mp4"))
"""

# ╔═╡ 39d9cd82-8dfb-48b0-827f-18dc187eda85
md"""
# Tracer distributions

$(LocalResource("../data/animations/cabbeling_600min/S_and_T_distributions.mp4"))
"""

# ╔═╡ e02371e2-64ce-4f90-bcaa-b0565366d14b
md"""
# Effective diffusivity

- Using the tracer percentile method (Sohail et al. (2020), Holmes et al. (2020)) we can estimate the effective vertical diffusivity as
```math
\kappa_{\mathrm{eff}} = \frac{\frac{\mathrm{d}}{\mathrm{d}t}\int_{V < V^{*}}\Theta(V, t) dV}{\frac{\partial \Theta}{\partial V}}.
```
- Figures below show that the effective diffusivity for the salinity tracer in the isothermal experiment are a good approximation to the parameterised value of ``\kappa_{S} = 1 \times 10^{-7}``.
$(LocalResource("iso_flux_diff.png"))
- The magnitude of the diffusivity in cabbeling case need further investigation.
$(LocalResource("flux_diff.png"))
"""

# ╔═╡ 27a2339e-300c-4e4e-a041-bc5251539c2a
md"""
# Energetics

- Following Winters et al. (1995) we can close the energy budget.
- In doing so we see that the conversion of APE to BPE is no longer irreversible in the cabbeling case.
$(LocalResource("irreversible_component.png"))
## Energy budget isothermal
$(LocalResource("iso_energies.png"))

## Energy budget cabbeling
$(LocalResource("cab_energies.png"))
"""

# ╔═╡ 0a303e91-42e1-4f44-99d7-38021f9f8ba6
md"""
# References
- A. Klocker et al. “Generation of the Internal Pycnocline in the Subpolar Southern Ocean by Wintertime Sea Ice Melting”. In: Journal of Geophysical Research: Oceans 128.3 (Mar. 2023). issn: 2169-9275. doi: 10.1029/2022JC019113.
- Jonas Nycander, Magnus Hieronymus, and Fabien Roquet. “The nonlinear equation of state of sea water and the global water mass distribution”. In: Geophysical Research Letters 42.18 (Sept. 2015), pp. 7714–7721. issn: 19448007. doi: 10.1002/2015GL065525.
- Fabien Roquet et al. “Unique thermal expansion properties of water key to the formation of sea ice on Earth”. In: Science Advances 8.46 (2022), p. 793. issn: 23752548. doi: 10.1126/ sciadv.abq0793.
- Winters, Kraig B., et al. "Available potential energy and mixing in density-stratified fluids." Journal of Fluid Mechanics 289 (1995): 115-128.
- Hughes, Graham O., Andrew Mc C. Hogg, and Ross W. Griffiths. "Available potential energy and irreversible mixing in the meridional overturning circulation." Journal of Physical Oceanography 39.12 (2009): 3130-3146.
"""

# ╔═╡ f941fdfb-d70a-473a-95ba-56f3aa53b477
# ╠═╡ disabled = true
#=╠═╡
md"""
# Effects of double diffusion
$(LocalResource("../data/animations/cabbeling_DD_600min_densityratio100/density.mp4"))
"""
  ╠═╡ =#

# ╔═╡ 25d49087-ab61-4b66-8981-21337448466f
# ╠═╡ disabled = true
#=╠═╡
md"""
$(LocalResource("../data/animations/cabbeling_DD_600min_densityratio100/tracers.mp4"))
"""
  ╠═╡ =#

# ╔═╡ bb1a2091-1087-4fb1-8904-145f8fd84fd2
# ╠═╡ disabled = true
#=╠═╡
md"""
$(LocalResource("../data/animations/cabbeling_DD_600min_densityratio100/hovmoller_∫ₐκᵥ.png"))
"""
  ╠═╡ =#

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
# ╟─39ba4eb8-4040-480a-9ea5-08d538f1f2ea
# ╟─ec83df92-8beb-4808-b1ac-8db61b81bae6
# ╟─e053dce3-1d3b-4012-afbe-65b1d754d8cb
# ╟─4c7769d2-852f-4c66-ad6c-31f453aedda9
# ╟─e26daadf-0134-463a-a182-2a63e5c73276
# ╟─0bffc2f9-1aa6-4cf6-9e3f-4010c97a7052
# ╟─3dc7a488-12a3-41e7-a4bf-a715072c800d
# ╟─ae8f6de0-7960-4beb-90af-7ea8e479cb56
# ╟─e46fcd1b-ccce-4089-8b49-66e93661fdc1
# ╟─f82ebb94-c94d-42cc-bd16-92b123ba6498
# ╟─bafe1643-5d0e-4d16-87ca-5edffef41c8c
# ╟─39d9cd82-8dfb-48b0-827f-18dc187eda85
# ╟─e02371e2-64ce-4f90-bcaa-b0565366d14b
# ╟─27a2339e-300c-4e4e-a041-bc5251539c2a
# ╟─0a303e91-42e1-4f44-99d7-38021f9f8ba6
# ╟─f941fdfb-d70a-473a-95ba-56f3aa53b477
# ╟─25d49087-ab61-4b66-8981-21337448466f
# ╟─bb1a2091-1087-4fb1-8904-145f8fd84fd2
# ╟─9feac837-c922-41b3-9e14-ba704aa2386b
