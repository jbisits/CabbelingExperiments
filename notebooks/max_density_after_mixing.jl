### A Pluto.jl notebook ###
# v0.19.29

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

# ╔═╡ 0f62e010-a6be-4e52-a10b-1773ce8941a9
@bind ulic Select(["stable", "cabbeling"])

# ╔═╡ 9950ba7b-5ea8-41bf-b721-088a52255ef9
S₀ᵘ = Dict("stable" => 34.551, "cabbeling" => 34.58)

# ╔═╡ 39d92993-9064-4549-9c25-48207432b052
begin
	S_star, Θ_star = 34.7, 0.5
	Θᵘ = -1.5
	slope = (Θᵘ - Θ_star) / (S₀ᵘ[ulic] - S_star)
	S_mix = range(S₀ᵘ[ulic], S_star, step = 0.000001)
	Θ_mix = @. Θᵘ - (slope) * (S₀ᵘ[ulic] - S_mix)
	ρ_mix = gsw_rho.(S_mix, Θ_mix, 0)
	max_rho, max_rho_idx = findmax(ρ_mix)
	Δρ_mix = max_rho - gsw_rho(S_star, Θ_star, 0)
	md"""
	Find that the gain in density is $Δρ_mix with the maximum density being $max_rho.
	This is at salinity $(round(S_mix[max_rho_idx], digits = 3))gkg⁻¹ and temperature $(round(Θ_mix[max_rho_idx], digits = 2))°C. So not at midpoint.
	"""
end 

# ╔═╡ 3f831fa4-d3a3-4059-aecf-1b693a337d99
let
	T₀ᵘ = -1.5
	S_star, Θ_star = 34.7, 0.5
	S_s, Θ_s = S₀ᵘ[ulic], T₀ᵘ
	Θ_c = T₀ᵘ
	N = 2000
	S_range, Θ_range = range(34.52, 34.72, length = N), range(-2, 1, length = N)
	S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
	ρ = gsw_rho.(S_grid, Θ_grid, 0)
	ρ_star = gsw_rho(S_star, Θ_star, 0)
	α_star = gsw_alpha(S_star, Θ_star, 0)
	β_star = gsw_beta(S_star, Θ_star, 0)
	ρ_s = gsw_rho(S_s, Θ_s, 0)
	find_Θ = findfirst(Θ_range .> -1.5)
	find_S = findfirst(ρ[find_Θ, :] .> ρ_star)
	S_iso, Θ_iso = S_range[find_S], Θ_range[find_Θ]
	gsw_rho(S_iso, Θ_iso, 0)
	αₗ, βₗ = gsw_alpha(S_star, Θ_star, 0), gsw_beta(S_star, Θ_star, 0)
	m = βₗ / αₗ
	Θ_linear = @. Θ_star + m * (S_range - S_star)
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1];
			  title = "Water masses in salinity-temperature space",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], color = :black, linewidth = 0.4, labelsize = 18, label = "Density at deep water")
	#lines!(ax, S_range, Θ_linear, color = :blue, linestyle = :dash)
	scatterlines!(ax, S_profile, T_profile, color = :purple, label = "Initial profile")
	scatter!(ax, [S_star], [Θ_star]; color = :red, label = "Lower layer")
	lines!(ax, S_mix, Θ_mix, color = :blue, label = "mixed water")
	scatter!(ax, S_mix[max_rho_idx], Θ_mix[max_rho_idx], color = :green, label = "maximum mixed density")
	contour!(ax, S_range, Θ_range, ρ', levels = [max_rho], color = :green, label = "maximum density")
	axislegend(ax, position = :lt)
	fig
	#save("initialTSprofiles.png", fig)
end

# ╔═╡ Cell order:
# ╟─f6dcc214-63e5-11ee-1a87-9ff13638a4b5
# ╟─99044830-2209-4f84-bade-b8f03f8af0ed
# ╟─0f62e010-a6be-4e52-a10b-1773ce8941a9
# ╟─9950ba7b-5ea8-41bf-b721-088a52255ef9
# ╟─39d92993-9064-4549-9c25-48207432b052
# ╠═3f831fa4-d3a3-4059-aecf-1b693a337d99
