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

# ╔═╡ 943cf012-755c-11ef-3493-99307e182e96
begin
	using Pkg
	Pkg.activate("..")
	using JLD2, CairoMakie, PlutoUI, StatsBase, Dates, GibbsSeaWater, Statistics	
end

# ╔═╡ b2ed38aa-523a-4edf-897a-5cf8c5f0ce33
md"""
# Cabbeling Energetics

This is a follow on from the energy budget analysis in `ha_diff_and_flux.jl`.
There I was looking at the energy budget and ways to isolate the effect of cabbeling on the BPE.
The main thing I ran into was *is it reasonable to claim we have reversible exhange between ``APE`` and ``BPE``?*
I know this is what the figures show but an exciting thing to consider is an alternated background state where mixing has been taken into account.

Recently I have been wondering if a Ramsdorf scheme might be able to give an approximation to the sorted then mixed state.
Questions arise from this instantly, one that I cannot figure out is *should the mixing take place on the adiabatically resorted profile?*
I suspect yes but how does this partition the ``BPE`` and does this mean the adiabatically resorted ``BPE`` is reversible and it is the sorted, mixed-sorted (because the point of mixing would be to resort in some way) that is only altered due to irreversible processes?
This seems reasonable to me but likely have to look at whether this makes sense mathematically.
"""

# ╔═╡ f0ebfa39-3959-4112-8575-e6238a73bfdc
begin

	bflux = load("buoyancy_flux.jld2")
	∫αΘw = bflux["∫αΘw"]
	∫βSw = bflux["∫βSw"]
	∫αΘw_interp = 0.5 * (∫αΘw[1:end-1] .+ ∫αΘw[2:end])
	∫βSw_interp = 0.5 * (∫βSw[1:end-1] .+ ∫βSw[2:end])
	bflux_keys = ["∫αΘw", "∫βSw"]
	
	bflux_face = load("buoyancy_flux_interp_face.jld2")
	time = bflux_face["time"]
	∫gρw = bflux_face["∫gρw"]
	∫Sw = bflux_face["∫Sw"]
	∫Θw = bflux_face["∫Θw"]
	βFₛ = bflux_face["βFₛ"]
	αFₜ = bflux_face["αFₜ"]
	energetics = load("cabbeling_energetics.jld2")
	pe = energetics["∫Ep_zref0"]
	bpe = energetics["∫Ebz✶_zref0"]
	ε = energetics["∫ε"]
	ek = energetics["∫Ek"]
	ρ₀ = energetics["ρ₀"]
	g = energetics["g"]
	z = energetics["z"]

	ape = pe .- bpe
	dₜek = diff(ek) ./ diff(time)
	dₜpe = diff(pe) ./ diff(time)
	dₜbpe = Φd = diff(bpe) ./ diff(time)
	dₜape = diff(ape) ./ diff(time)
	time_interp = 0.5 * (time[1:end-1] .+ time[2:end])
	∫gρw_interp = Φz = 0.5 * (∫gρw[1:end-1] .+ ∫gρw[2:end]) ./ ρ₀
	∫Sw_interp = 0.5 * (∫Sw[1:end-1] .+ ∫Sw[2:end])
	∫Θw_interp = 0.5 * (∫Θw[1:end-1] .+ ∫Θw[2:end])
	ε_interp = 0.5 * (ε[1:end-1] .+ ε[2:end])
	Φi = dₜpe .- Φz
	ΔV = load("cabbeling_fluxes_and_diff_longer_run.jld2", "ΔV")
	βFₛ = 0.01 * βFₛ / (g * ΔV) # Scale appropriately --- need to check!
	αFₜ = 0.01 * αFₜ / (g * ΔV)
	md"""

	The data is currently saved across a few files so here is a breakdown of what is in which.
	Where appropriate output has been interpolated in time.
	Data loaded from `buoyancy_flux_interp_face.jld2`:
	
	$([string(k)*", " for k ∈ keys(bflux_face)])

	data loaded from `buoyancy_flux.jld2`:

	$([string(k)*", " for k ∈ keys(bflux)])

	data loaded from `cabbeling_energetics.jld2`:

	$([string(k)*", " for k ∈ keys(energetics)]).
	"""
end

# ╔═╡ 5bf45be1-f79d-48e2-a472-f6f1b8935bf3
md"""
## Energy budget

This has already been shown in `ha_diff_and_flux.jl` but just for thoroughness, in our simulations (closed domain) we should have time change in kinetic energy be due to dissipation and buoyancy flux,
```math
\frac{\mathrm{d}}{\mathrm{d}t} E_{k} = - ε - g∫ρw\mathrm{d}V.
```
"""

# ╔═╡ 38a60f55-3756-4554-9e09-08680e6cb8f5
@bind ek_window PlutoUI.Slider(2:length(time_interp), default = 200)

# ╔═╡ 2ad9595c-62e1-41e1-8905-20e208150629
let
	rhs = -ε_interp .- ∫gρw_interp
	fig, ax = lines(time_interp[1:ek_window], dₜek[1:ek_window], label = "dₜEk")
	lines!(time_interp[1:ek_window], rhs[1:ek_window], label = "-ε - g∫ρwdV")
	MSE = mean((dₜek[1:ek_window] .- rhs[1:ek_window]).^2)
	ax.title = "Volume intergated kinetic energy change"
	ax.subtitle = "MSE = $(MSE)"
	ax.ylabel = "Watts"
	axislegend(ax, position = :rt)
	ax2 = Axis(fig[2, 1], xlabel = "time (s)", ylabel = "Watts", title = "Absolute error")
	abs_err = abs.(dₜek[1:ek_window] .- rhs[1:ek_window])
	lines!(ax2, time_interp[1:ek_window], abs_err)
	fig
end

# ╔═╡ 2ce30c05-8923-49df-ba80-945d85636507
md"""
## Total potential, background potential and available potential energy

The calculations of the total potential and background potential energy have been reference to ``z = 0``.
This means in the `sort`ing of the ``BPE``, the lightest elements are at the top (``z = 1``) and the heaviest elements are at the bottom (``z = 0``).
For the ``PE`` this meant shifting the `z` grid so ``z ∈ [0, 1]`` rather than ``[-1, 0]``.
The ``APE`` is still the difference between the two.

The only thing I do not understand is why initially ``BPE > PE``?
There is a small amount of noise at the interface but sorting the density profile still looks different and, importantly, more stable.
What is going on there?
Otherwise this looks better referenced from ``z = 0``.
"""

# ╔═╡ e61e43f7-d9b8-4a8d-8a3d-68bdf41c517f
@bind pe_window PlutoUI.Slider(2:length(time_interp), default=200)

# ╔═╡ c26744ec-7f38-4ab3-b68f-c99ab313a53e
let
	fig, ax = lines(time[1:pe_window], pe[1:pe_window], label = "PE")
	lines!(ax, time[1:pe_window], bpe[1:pe_window], label = "BPE")
	ax.title = "Potential and background potential energies"
	ax.ylabel = "Joules"
	axislegend(ax, position = :rb)
	ax2 = Axis(fig[2, 1], xlabel = "time", ylabel = "Joules", title = "Available potential energy")
	lines!(time[1:pe_window], ape[1:pe_window], color = :red, label = "APE")
	axislegend(ax2, position = :rt)
	fig
end

# ╔═╡ 8d4756a8-06a4-4436-9a27-e1702a65bc14
let
	fig, ax = lines(time_interp[1:pe_window], dₜpe[1:pe_window], label = "dₜPE")
	lines!(ax, time_interp[1:pe_window], dₜbpe[1:pe_window], label = "dₜBPE")
	ax.title = "Time derivatie of potential and background potential energies"
	ax.ylabel = "Watts"
	axislegend(ax, position = :rt)
	ax2 = Axis(fig[2, 1], xlabel = "time", ylabel = "Watts", title = "Time derivative of available potential energy")
	lines!(ax2, time_interp[1:pe_window], dₜape[1:pe_window], label = "dₜAPE", color = :red)
	axislegend(ax2, position = :rt)
	fig
end

# ╔═╡ 786f2208-39ff-4fd2-8363-b6aa41c42e7b
let
	fig, ax = lines(time_interp[1:pe_window], Φz[1:pe_window], label = "Φz")
	lines!(ax, time_interp[1:pe_window], Φd[1:pe_window], label = "Φd")
	lines!(ax, time_interp[1:pe_window], Φi[1:pe_window], label = "Φi")
	lines!(ax, time_interp[1:pe_window], ε_interp[1:pe_window], label = "ε")
	ax.title = "Enegry exchanges"
	ax.xlabel = "time (s)"
	ax.ylabel = "Watts"
	axislegend(ax, position = :rt)
	fig
end

# ╔═╡ c3585182-3b5f-4a5d-8616-5c06f35c8811
md"""
Computing ``Φ_{i} = \mathrm{d}_{t}PE - Φz`` as a residual seems questionable here.
"""

# ╔═╡ b0efe882-59a1-4319-99f2-04fb59f062c4
md"""
## Salinity and temperature fluxes

The goal here is to find the fluxes due the diffusive fluxes of salinity and temperature, use these as a proxy for density flux that should be conserved to isolate the non-linear contribution.
"""

# ╔═╡ c73bdae3-754c-405c-8599-7f0767fed5c9
let
	fig, ax = lines(time_interp, βFₛ, label = "βFₛ", color = :blue)
	ax.xlabel = "time"
	axislegend(ax, position = :rb)
	ax2 = Axis(fig[1, 2], xlabel = "time (s)")
	lines!(ax2, time_interp, αFₜ, label = "αFₜ", color = :red)
	axislegend(ax2)
	ax3 = Axis(fig[2, :])
	S_and_T_fluxes = (αFₜ .- βFₛ)
	lines!(ax3, time_interp, S_and_T_fluxes, label = "∫αΘw - ∫βSw", color = :purple)
	lines!(ax3, time_interp, dₜbpe, label = "dₜbpe", color = :orange)
	axislegend(ax3)
	hidexdecorations!(ax3, grid = false, ticks = false)
	ax4 = Axis(fig[3, :], xlabel = "time (s)")
	lines!(ax4, time_interp,  dₜbpe .- S_and_T_fluxes, color = :green, label = "dₜbpe - (∫αΘw - ∫βSw)", title = "Cabbeling contribution (?)")
	axislegend(ax4)
	fig
end

# ╔═╡ ffd3a976-1234-4758-8ec2-643be78af3f7
TableOfContents(title = "Energetics analysis")

# ╔═╡ Cell order:
# ╟─943cf012-755c-11ef-3493-99307e182e96
# ╟─b2ed38aa-523a-4edf-897a-5cf8c5f0ce33
# ╠═f0ebfa39-3959-4112-8575-e6238a73bfdc
# ╟─5bf45be1-f79d-48e2-a472-f6f1b8935bf3
# ╟─38a60f55-3756-4554-9e09-08680e6cb8f5
# ╟─2ad9595c-62e1-41e1-8905-20e208150629
# ╟─2ce30c05-8923-49df-ba80-945d85636507
# ╟─e61e43f7-d9b8-4a8d-8a3d-68bdf41c517f
# ╟─c26744ec-7f38-4ab3-b68f-c99ab313a53e
# ╟─8d4756a8-06a4-4436-9a27-e1702a65bc14
# ╟─786f2208-39ff-4fd2-8363-b6aa41c42e7b
# ╟─c3585182-3b5f-4a5d-8616-5c06f35c8811
# ╟─b0efe882-59a1-4319-99f2-04fb59f062c4
# ╟─c73bdae3-754c-405c-8599-7f0767fed5c9
# ╟─ffd3a976-1234-4758-8ec2-643be78af3f7
