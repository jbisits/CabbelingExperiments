### A Pluto.jl notebook ###
# v0.19.35

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

# ╔═╡ 99044830-2209-4f84-bade-b8f03f8af0ed
begin
	using Pkg
	Pkg.activate("..")
	using CairoMakie, GibbsSeaWater, PlutoUI
end

# ╔═╡ f6dcc214-63e5-11ee-1a87-9ff13638a4b5
md"""
# Maximum density after mixing

In `initial_stepchange.jl` there is a figure I wish to explore more.
When computing the maximum density from the mixed water this also coincided exactly with the **linearised density at the post mixing maximum density.**
Is this luck or something that could be further explored?
"""

# ╔═╡ 3e1f2eb1-45cd-44bc-8fec-85acd0abda4d
begin
	shallow_salinity = @bind S₀ᵘ PlutoUI.Slider(range(34.551, 34.6, step = 0.001))
	nothing
end

# ╔═╡ 39d92993-9064-4549-9c25-48207432b052
begin
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
			  title = "Mixing water masses to find new maximum density",
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
	Fix the deep water at ``S = `` $(S_star)gkg⁻¹ ``\Theta = `` $(Θ_star)°C.
	We then set the shallow water at ``\Theta = `` $(Θᵘ)°C and vary ``S`` between stable to cabbeling and isohaline with the slider

	``S`` = $(shallow_salinity) = $(S₀ᵘ)gkg⁻¹.

	The maximum density after mixing is $(round(max_rho, digits = 5)) kgm⁻³ which is a gain of $(round(Δρ_mix, digits = 5))kgm⁻³.
	The new maximum density is at salinity $(round(S_mix[max_rho_idx], digits = 3))gkg⁻¹ and temperature $(round(Θ_mix[max_rho_idx], digits = 2))°C.

	Slope of mixing line is $(m_initial) and linearised density slope is $(m).

	As we get closer to the initial shallow water being at the same density as the deep water, salinity and temperature approach the midpoint values.

	$(fig)
	"""

end

# ╔═╡ Cell order:
# ╟─f6dcc214-63e5-11ee-1a87-9ff13638a4b5
# ╟─99044830-2209-4f84-bade-b8f03f8af0ed
# ╟─3e1f2eb1-45cd-44bc-8fec-85acd0abda4d
# ╟─39d92993-9064-4549-9c25-48207432b052
