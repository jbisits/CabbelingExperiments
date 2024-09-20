### A Pluto.jl notebook ###
# v0.19.46

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

# ╔═╡ 7aee38b0-762d-11ef-06a5-3b50328ccb52
begin
	using Pkg
	Pkg.activate("..")
	using GibbsSeaWater, PlutoUI, CairoMakie, StatsBase
end

# ╔═╡ 780b30dc-cd32-4867-86c6-d8aefff90a16
md"""
# Estimating the lowest background potential energy state

The results of project two show that it is not sufficient to adiabatically sort Boussinesq fluid elements from a grid to get a state of minimal potential energy for a system.
It will be when the linear equation of state assumption is met but introducing non-linearity means density is not conserved when mixed so we have
```math
\frac{∂ρ}{∂t} + \mathbf{u} ⋅ ∇ρ = κ∇^2ρ + F_{NL},
```
or something like it, in place of equation (6) in Winters et al. (1995) (the ``F_{NL}`` term is supposed to be a positive definite source term of density due to non-lienarity).
If the above is permissable, then the evolution of the ``BPE`` following Winters (1995) is
```math
\frac{\mathrm{d}}{\mathrm{d}t}E_{b} = g∫_{V}\left(-\mathbf{u} ⋅ ∇ρ + κ∇^2ρ + F_{NL}\right)z✶\mathrm{d}V.
```
Linearity allows us to break this up so we end up with the result of Winters (1995), equation (18) or a previous variant of it, and the non-linear term
```math
\frac{\mathrm{d}}{\mathrm{d}t}E_{b} = κg∫_{V}-\frac{\mathrm{d}z✶}{\mathrm{d}ρ}|∇ρ|^2\mathrm{d}V  + g∫_{V}F_{NL} z✶ \mathrm{d}V.
```

More could/should be done mathematically to show something about this.

What we would like to be able to do in this current project is to say something about the term ``g∫_{V}F_{NL} z✶ \mathrm{d}`` from the residual of ``BPE`` and the *linear* part (currently I am working on this using the fluxes of salt and temperature).

Maybe there is also something that could be done in the estimating of the ``BPE`` though.
With my simulations I know what the minimal PE state --- it is when all deep water is brought to the maximum density and it is only the background diffusivity that alters this profile.
In essence I think *this* is the background state, one where the process is irreversible.

One thing I thought of is a Rahmsdorf scheme *on the adiabatically sorted profile*.
This would need to be done in a conservative manner but we could look at a sorted density profile, get the ``S`` and ``Θ`` permutation from the sorted density and apply the Rahmsdorf scheme:
- look at two vertically stacked levels and mix the ``S`` and ``\Theta`` conservatively;
- check ``Δρ`` between these mixed levels and the next level;
- if unstable mix again else leave there.

I am not sure how this will work but interested to try.
"""

# ╔═╡ ff1abe25-6874-4968-a707-0361a0a9be56
md"""
## Synthetic profile

I will use the initial condition from my cabbeling and isothermal experiments.
At this stage no evolution -- I just want to see what happens to the estimated background state.
Note I reference everything here from ``z = 0``.
With the noise I get my `sort`ed state looking slightly different which helps.
"""

# ╔═╡ 790c24ee-702b-42f5-b6aa-cd9d4c461002
@bind expt Select(["cabbeling", "isothermal"])

# ╔═╡ 6ab10d37-1710-4060-8a7c-671ba425aac9
function Rahmsdorf(S_profile, Θ_profile, σ_profile, σ_perm; pᵣ = 0)

	mixed_sorted_σ_profile = similar(S_profile)
	sorted_σ_profile = σ_profile[σ_perm]
	sorted_S = S_profile[σ_perm]
	sorted_Θ = Θ_profile[σ_perm]
	
	for i ∈ length(S_profile):-1:3

		Sₘ = mean(sorted_S[i:i-1])
		Θₘ = mean(sorted_Θ[i:i-1])
		σₘ = gsw_rho.(Sₘ, Θₘ, pᵣ)

		if σₘ > sorted_σ_profile[i-2]

			j = 1
			unstable = true
			while unstable
				Sₘ_ = mean(sorted_S[i:i-j])
				Θₘ_ = mean(sorted_Θ[i:i-j])
				σₘ_ = gsw_rho.(Sₘ_, Θₘ_, pᵣ)
				unstable = σₘ_ > sorted_σ_profile[i-j-1] && j > 2 ? true : false
				mixed_sorted_σ_profile[i-j-1] = σₘ_
				j += 1
			end
			
		else
			mixed_sorted_σ_profile[i-2] = sorted_σ_profile[i-2]
		end
	end
	mixed_sorted_σ_profile[end-1:end] = sorted_σ_profile[end-1:end]
	return mixed_sorted_σ_profile
end

# ╔═╡ 3de183d1-d3a7-410e-986d-c5beaf592be9
gsw_rho.(0.5*(34.58 + 34.7), 0.5*(-1.5 + 0.5), 0) > gsw_rho(34.7, 0.5, 0)

# ╔═╡ 113788b0-6adc-4d27-a19f-af95755b27ea
begin
	Sₗ, Θₗ = 34.7, 0.5
	Θᵤ = expt == "cabbeling" ? -1.5 : 0.5
	Sᵤ = expt == "cabbeling" ? 34.58 : 34.69431424
	pᵣ = 0
	
	N = 1000
	z = range(0, 1, N)
	interface_location = 0.5
	find_interface = findall(interface_location - 0.02 .< z .< interface_location + 0.02)
	S_profile = vcat(fill(Sₗ, Int(N / 2)), fill(Sᵤ, Int(N / 2)))
	S_profile[find_interface] .+= 2e-4 * randn(length(find_interface))
	Θ_profile = vcat(fill(Θₗ, Int(N / 2)), fill(Θᵤ, Int(N / 2)))

	σ_profile = gsw_rho.(S_profile, Θ_profile, pᵣ)
	sorted_σ_perm = sortperm(σ_profile, rev = true)
	sorted_σ_profile = σ_profile[sorted_σ_perm]
end

# ╔═╡ c85f96bd-6ff2-4a97-adb4-7f49b6daa6fd
test = Rahmsdorf(S_profile, Θ_profile, σ_profile, sorted_σ_perm)

# ╔═╡ a1cd1c2b-207e-4a90-aa07-a25b3f926b2a
let
	fig, ax = lines(σ_profile, z, label = "In-situ")
	ax.title = "Initial density profile"
	ax.xlabel = "σ"
	ax.ylabel = "z"
	lines!(ax, sorted_σ_profile, z, label = "Sorted")
	lines!(ax, test, z)
	axislegend(ax, position = :rt)
	fig
end

# ╔═╡ 1769fdbb-50b4-45e3-9450-221ce20f254e
TableOfContents()

# ╔═╡ Cell order:
# ╟─7aee38b0-762d-11ef-06a5-3b50328ccb52
# ╟─780b30dc-cd32-4867-86c6-d8aefff90a16
# ╟─ff1abe25-6874-4968-a707-0361a0a9be56
# ╟─790c24ee-702b-42f5-b6aa-cd9d4c461002
# ╟─6ab10d37-1710-4060-8a7c-671ba425aac9
# ╟─3de183d1-d3a7-410e-986d-c5beaf592be9
# ╟─113788b0-6adc-4d27-a19f-af95755b27ea
# ╟─a1cd1c2b-207e-4a90-aa07-a25b3f926b2a
# ╟─c85f96bd-6ff2-4a97-adb4-7f49b6daa6fd
# ╟─1769fdbb-50b4-45e3-9450-221ce20f254e
