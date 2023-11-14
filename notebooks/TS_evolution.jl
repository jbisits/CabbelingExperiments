### A Pluto.jl notebook ###
# v0.19.32

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

# ╔═╡ ced5de2d-6e19-4e32-832c-1dcdcca82642
begin
	using Pkg
	Pkg.activate("..")
	using TwoLayerDirectNumericalShenanigans
	using Oceananigans.Fields
	using PlutoUI, GibbsSeaWater
	using CairoMakie
	using SpecialFunctions: erf
end

# ╔═╡ 3d90f332-f3f1-45f1-bfea-3fa4d1c648af
begin
	import TwoLayerDirectNumericalShenanigans.erf_tracer_solution
	erf_tracer_solution(z, Cₗ::Number, ΔC::Number, κ::Number, time, interface_location) =
    Cₗ + 0.5 * ΔC * (1 + erf((z - interface_location) / sqrt(4 * κ * time)))
end

# ╔═╡ faa94bce-2a04-11ee-39f3-518323d8ad0f
md"""
# Analytic evolution of temperature and salinity from step function

A 1D profile of temperature (or salinity) initially as a step function will evolve according to
```math
	\begin{aligned}
	\frac{\partial \Theta}{\partial t} &= \kappa_{\Theta}\frac{\partial^2 \Theta}{\partial z^2} \\
	\Theta(z, 0) &= \begin{cases} 0 & z < 0 \\ 1 & z > 0. \end{cases}
	\end{aligned}
```

This can be solved analytically to give the evolution of the tracer in time
```math
	\Theta(z, t) = \frac{1}{2}\left(1 + \mathrm{erf}\left(\frac{z}{\sqrt{4\kappa_{\Theta}t}}\right)\right).
```

The location of the step change and the magnitude of difference can be added so that the initial condition becomes
```math
	\Theta(z, 0) = \begin{cases} \Theta^{*} & z < -a \\ \Theta^{*} + \Delta\Theta & z > -a \end{cases}
```
which is the Heaviside function ``\Theta^{*} + \Delta\Theta H\left(z + a\right).``
This is the kind of step change profile we have looked already and it can be modelled analytically according to
```math
	\Theta(z, t) = \Theta^{*} + \frac{\Delta\Theta}{2}\left(1 + \mathrm{erf}\left(\frac{z + a}{\sqrt{4\kappa_{\Theta}t}}\right)\right)
```
with a similar expression for salinity.
So we set the lower layer temperature and salinity as ``\left(S^{*}, \Theta^{*}\right)`` then pick an interface and ``\left(\Delta S, \Delta\Theta\right)`` that the upper layer is and we can evolve salinity and temperature and see how density evolves due to mixing and different initial conditions.

We aim to do something, maybe even use these analytic models, in the DNS.
Along with looking at the DNS, we will use a 1D model to investigate how convective adjustment schemes compare to the DNS.
Different salinity and temperature diffusivities can also be explored using these models as below.

## 1D model

The 1D model has ranges over ``z \in \left(-500m, 0\right)``, has the step change interface at ``z = -100m``, and has a vertical resolution of ``1m.``
The only dynamics come from the convective adjustment scheme.

In all plots below I set the lower layer as ``\left(S^{*}, \Theta^{*}\right) = \left(34.7, 0.5\right)`` and use a stable and cabbeling initial condition to explore what is happening.

### Equal diffusivity for salinity and temperature

Running the model with equal diffusivity for temperature and salinity triggered the convective adjustment scheme at the **nine minute mark** in the simulation.
From this analytic solution the first bulge in the density appears (as in appears visually) around the 5 minute mark.

"""

# ╔═╡ 894a827f-ba32-4e79-89ce-d9663288293b
@bind upper_layer Select(["Stable", "Cabbeling"])

# ╔═╡ da22a77d-d1e0-4da6-832b-8c1134d774b4
begin
	S_star, Θ_star = 34.7, 0.5
	ΔΘ = -2
	S_stable, S_cab = 34.551, 34.58
	ΔSₛ = S_stable - S_star
	ΔS_c= S_cab - S_star
	ΔS = upper_layer == "Stable" ? ΔSₛ : ΔS_c
	nothing
end

# ╔═╡ 8cf080c8-70bf-48f7-9cb4-dbdb4a492941
begin
	time = @bind times PlutoUI.Slider(0.00001:10000.00001)
	κₛ = @bind κₛ PlutoUI.Slider([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9], show_value = true)
	κₜ = @bind κₜ PlutoUI.Slider([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9], show_value = true)
	md"""
	### Diffusivities
	Salt = $(κₛ) m²s⁻¹

	Temperature = $(κₜ) m²s⁻¹

	### Time
	t = $(time)
	"""
end

# ╔═╡ 83147607-4047-4bd6-8eb6-a8e17fdd0c84
@bind zoom Select(["Full profile", "Zoom on interface"])

# ╔═╡ 0194e952-56a2-4718-8ac2-5934f044cce3
let

	z = range(-500, 0, length = 500)
    interface_location = -100.0
	S = erf_tracer_solution.(z, S_star, ΔS, κₛ, times, interface_location)
	T = erf_tracer_solution.(z, Θ_star, ΔΘ, κₜ, times, interface_location)
	Sᵢ = erf_tracer_solution.(z, S_star, ΔS, κₛ, 1e-5, interface_location)
	Tᵢ = erf_tracer_solution.(z, Θ_star, ΔΘ, κₜ, 1e-5, interface_location)
	σ₀ᵢ = gsw_rho.(Sᵢ, Tᵢ, 0)
	σ_limits = (σ₀ᵢ[end] - 0.004, σ₀ᵢ[1] + 0.04)
	σ₀ = gsw_rho.(S, T, 0)

	fontsize = 22
	labelsize = 16
	fig = Figure(size = (900, 400); fontsize)
	ax = [Axis(fig[1, i]) for i ∈ 1:2]
	lines!(ax[1], S, z; color = (:blue, 0.5), label = "Salinity")
	xlims!(ax[1], 34.52, 34.72)
	ax[1].xlabel = "Salinity (gkg⁻¹)"
	ax[1].xlabelcolor = :blue
	ax[1].xticklabelcolor = :blue
	ax[1].ylabel = "z (m)"
	ax2 = Axis(fig[1, 1];
	           xaxisposition = :top,
	           xticklabelcolor = :red,
	           xlabel = "Θ (°C)",
	           xlabelcolor = :red,
	           title = "Temperature and salinity profiles\nat t = $(round(times/60; digits = 2)) mins")
	lines!(ax2, T, z; color = (:red, 0.5), label = "Temeperature")
	lines!(ax[2], σ₀, z; color = (:black, 0.5))
	ax[2].title = "Density at t = $(round(times/60; digits = 2)) mins."
	ax[2].xlabel = "σ₀ (kgm⁻³)"
	ax[2].xticklabelrotation = π/4
	xlims!(ax[2], σ_limits)
	axislegend(ax2; labelsize)
	axislegend(ax[1], position = :lb; labelsize)
	hideydecorations!(ax[2], grid = false)

	if isequal(zoom, "Zoom on interface")
		ylims!(ax[1], -150, -50)
		ylims!(ax[2], -150, -50)
		ylims!(ax2, -150, -50)
	end
	linkyaxes!(ax[1], ax[2])
	colsize!(fig.layout, 1, Relative(3/5))
	fig

end

# ╔═╡ 56fb1d65-2bcf-45ca-88f8-f310c9d992b4
md"""
### Different diffusivity for salinity and temperature

With different diffusivity for salinity and temperature (``\kappa_{S}`` two orders of magnitude smaller as is thought to be the case in Southern Ocean) the first bulge in density is around 3 minute mark.
Currently I do not have the model running with different salinity and temperture tracers but this will be an interesting thing to look at - does the slower salt diffusivity allow the temperature to more quickly mix and create an instability?

As well as this the "stable" initial condition sees the same instability form as in the "cabbeling" initial condition.
It must be that when the salt is not able to be difffused at the same rate to maintain the stratification, the temperature diffusion creates the cabbeling instability.
Contrast this with the situation above with equal salinity and temperature diffusivity where no instability forms (have not tested to infinity..).
Interesting!

As this analytic solution is stepped forward, we can also see decreas in density at the interface closer to the upper layer.
I am not sure why this is the case but will be somehting to look at more, to see how this affects the convective adjustment scheme, once the 1D model can run with different ``S`` and ``\Theta`` diffusivities.
"""

# ╔═╡ c8b85349-2358-43fc-8aa6-1e70adf4feb5
md"""
## DNS

The DNS has an overall depth of ``1m`` and is ``0.1m`` in the horizontal directions.
The resolution is very high but not uniform over the domain as there is stretching at the bottom to reduce the computational cost.

At this stage I have not run any DNS with this configuration, so I have nothing to compare this to.
What I hope is to use the analytic model to create initial temperature and salinity profiles though this does has the limitaion that at the interface there will be a density perturbation.
This may mean that instead of this we set stable density then have a salinity pertubration above the interface to kick off the mixing - still thinking about what might be most appropriate.

### Equal diffusivity for salinity and temperature
"""

# ╔═╡ 828c6a69-a61f-4a98-ba06-e09661d8a973
begin
	timeDNS = @bind timeDNS PlutoUI.Slider(0.00001:100000.00001)
	κₛDNS = @bind κₛDNS PlutoUI.Slider([1e-5, 1e-6, 1e-7, 1e-8, 1e-9], show_value = true)
	κₜDNS = @bind κₜDNS PlutoUI.Slider([1e-5, 1e-6, 1e-7, 1e-8, 1e-9], show_value = true)
	md"""
	### Diffusivities
	Salt = $(κₛDNS) m²s⁻¹

	Temperature = $(κₜDNS) m²s⁻¹

	### Time
	t = $timeDNS
	"""
end

# ╔═╡ 7b5ff0ab-1fbd-4252-9d3a-27285c9b5dab
let

	z = range(-1, 0, length = 500)
    interface_location = -0.375
	S = TwoLayerDirectNumericalShenanigans.erf_tracer_solution.(z, S_star, ΔS, κₛDNS, timeDNS, interface_location)
	T = TwoLayerDirectNumericalShenanigans.erf_tracer_solution.(z, Θ_star, ΔΘ, κₜDNS, timeDNS, interface_location)
	σ₀ = gsw_rho.(S, T, 0)

	fontsize = 22
	labelsize = 16
	fig = Figure(size = (900, 400); fontsize)
	ax = [Axis(fig[1, i]) for i ∈ 1:2]
	lines!(ax[1], S, z; color = (:blue, 0.5), label = "Salinity")
	xlims!(ax[1], 34.52, 34.72)
	ax[1].xlabel = "Salinity (gkg⁻¹)"
	ax[1].xlabelcolor = :blue
	ax[1].xticklabelcolor = :blue
	ax[1].ylabel = "z (m)"
	ax2 = Axis(fig[1, 1];
	           xaxisposition = :top,
	           xticklabelcolor = :red,
	           xlabel = "Θ (°C)",
	           xlabelcolor = :red,
	           title = "Temperature and salinity profiles\nat t = $(round(timeDNS/60; digits = 4)) mins")
	lines!(ax2, T, z; color = (:red, 0.5), label = "Temeperature")
	lines!(ax[2], σ₀, z; color = (:black, 0.5))
	ax[2].title = "Density at t = $(round(timeDNS/60; digits = 4)) mins."
	ax[2].xlabel = "σ₀ (kgm⁻³)"
	ax[2].xticklabelrotation = π/4
	axislegend(ax2; labelsize)
	axislegend(ax[1], position = :lb; labelsize)
	hideydecorations!(ax[2], grid = false)

	linkyaxes!(ax[1], ax[2])
	colsize!(fig.layout, 1, Relative(3/5))
	fig

	fig2 = Figure(size = (500, 500); fontsize)
	ax2 = Axis(fig2[1, 1], title = "S-T evolution", xlabel = "Absolute salinity (gkg⁻¹)", ylabel = "Conservative temperature (°C)")
	N = 2000
	S_range, Θ_range = range(minimum(S), maximum(S), length = N), range(minimum(T), maximum(T), length = N)
	S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
	ρ = gsw_rho.(S_grid, Θ_grid, 0)
	ρ_star = gsw_rho(S_star, Θ_star, 0)
	contour!(ax2, S_range, Θ_range, ρ'; levels = [ρ_star], color = :red, label = "Isopycnal at deep layer")

	scatter!(ax2, S, T, σ₀, markersize = 6)

	md"""
	$(fig)
	$(fig2)
	"""

end

# ╔═╡ 70b786cc-81d9-4a25-a911-7adb78e8ae82
md"""
### Different diffusivity for salinity and temperature
"""

# ╔═╡ adaf329a-9316-4e0c-b7e3-cba513cebf68
md"""
# Numerical evolution of temperature and salinity
"""

# ╔═╡ f73d2ef6-fbca-4306-a838-0079ee34d6dd
begin
	output = joinpath(@__DIR__, "../1Dmodel/OneDModelOutput_cabbeling.jld2")
	S_ts, T_ts = FieldTimeSeries(output, "S"), FieldTimeSeries(output, "T")
	model_times_idx = @bind mt_idx PlutoUI.Slider(eachindex(S_ts.times))
	model_times = S_ts.times
	z_model = znodes(S_ts.grid, Center(), Center(), Center())
	nothing
end

# ╔═╡ e04e5e49-8b5b-4c01-9184-c371c3651cdd
md"""
t = $model_times_idx $(model_times[mt_idx] / 60) mins
"""

# ╔═╡ 5a020f83-0c13-4ca0-94fc-756c5b0494ba
@bind zoom_model Select(["Full profile", "Zoom on interface"])

# ╔═╡ 03494eb9-b3fd-4380-9492-014d4954b73a
let

	S = interior(S_ts, 1, 1, :, mt_idx)
	T = interior(T_ts, 1, 1, :, mt_idx)
	σ₀ = gsw_rho.(S, T, 0)

	fontsize = 22
	labelsize = 16
	fig = Figure(size = (900, 400); fontsize)
	ax = [Axis(fig[1, i]) for i ∈ 1:2]
	lines!(ax[1], S, z_model; color = (:blue, 0.5), label = "Salinity")
	xlims!(ax[1], 34.52, 34.76)
	ax[1].xlabel = "Salinity (gkg⁻¹)"
	ax[1].xlabelcolor = :blue
	ax[1].xticklabelcolor = :blue
	ax[1].ylabel = "z (m)"
	ax2 = Axis(fig[1, 1];
	           xaxisposition = :top,
	           xticklabelcolor = :red,
	           xlabel = "Θ (°C)",
	           xlabelcolor = :red,
	           title = "Temperature and salinity profiles\nat t = $(model_times[mt_idx] / 60) mins")
	lines!(ax2, T, z_model; color = (:red, 0.5), label = "Temeperature")
	lines!(ax[2], σ₀, z_model; color = (:black, 0.5))
	ax[2].title = "Density at t = $(model_times[mt_idx] / 60) mins."
	ax[2].xlabel = "σ₀ (kgm⁻³)"
	ax[2].xticklabelrotation = π/4
	#xlims!(ax[2], σ_limits)
	axislegend(ax2; labelsize)
	axislegend(ax[1], position = :lb; labelsize)
	hideydecorations!(ax[2], grid = false)

	if isequal(zoom_model, "Zoom on interface")
		ylims!(ax[1], -150, -50)
		ylims!(ax[2], -150, -50)
		ylims!(ax2, -150, -50)
	end
	linkyaxes!(ax[1], ax[2])
	colsize!(fig.layout, 1, Relative(3/5))
	fig
end

# ╔═╡ 090d2949-cf5a-41e4-8216-460d207fc0f5
TableOfContents(title = "Analytic and numeric temperature and salinity evolution")

# ╔═╡ Cell order:
# ╟─ced5de2d-6e19-4e32-832c-1dcdcca82642
# ╟─3d90f332-f3f1-45f1-bfea-3fa4d1c648af
# ╟─faa94bce-2a04-11ee-39f3-518323d8ad0f
# ╟─894a827f-ba32-4e79-89ce-d9663288293b
# ╟─da22a77d-d1e0-4da6-832b-8c1134d774b4
# ╟─8cf080c8-70bf-48f7-9cb4-dbdb4a492941
# ╟─83147607-4047-4bd6-8eb6-a8e17fdd0c84
# ╟─0194e952-56a2-4718-8ac2-5934f044cce3
# ╟─56fb1d65-2bcf-45ca-88f8-f310c9d992b4
# ╟─c8b85349-2358-43fc-8aa6-1e70adf4feb5
# ╟─828c6a69-a61f-4a98-ba06-e09661d8a973
# ╟─7b5ff0ab-1fbd-4252-9d3a-27285c9b5dab
# ╟─70b786cc-81d9-4a25-a911-7adb78e8ae82
# ╟─adaf329a-9316-4e0c-b7e3-cba513cebf68
# ╟─f73d2ef6-fbca-4306-a838-0079ee34d6dd
# ╟─e04e5e49-8b5b-4c01-9184-c371c3651cdd
# ╟─5a020f83-0c13-4ca0-94fc-756c5b0494ba
# ╟─03494eb9-b3fd-4380-9492-014d4954b73a
# ╟─090d2949-cf5a-41e4-8216-460d207fc0f5
