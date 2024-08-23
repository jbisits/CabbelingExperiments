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

# ╔═╡ b50f0e02-5ea4-11ef-0758-95c6ef87071e
begin
	using Pkg
	Pkg.activate("..")
	using JLD2, CairoMakie, PlutoUI, StatsBase, TimeSeries, Dates, GibbsSeaWater
end

# ╔═╡ d11f7aad-5cce-4957-a9cb-8366e2b219f1
md"""
# Horizontally averaged effective diffusivity and fluxes

One way to get rid of the spike in the effective diffusivity is to first horizontally average.
This gives behaviour like we would expect and perhaps taking a moving mean we can generate a good estimate using this diagnostic.

The effective diffusivity at each depth level of the *horizontally averaged salinity profile* is
```math
\kappa_{\mathrm{eff}}(z^{*}, t) = \frac{∂_{t}\int_{z < z^{*}} S \mathrm{d}z}{\partial S / \partial z}.
```
Integrating over all depths and dividing by total depth, which is 1m, to get a an average effective diffusivity for a *horizontally averaged profile*.
**Note:** the horizontally averaged profile is still **sorted**.

Tracer fluxes from horizontally averaged, sorted profiles are
```math
F_{C} = ∂_{t}\int_{z < z^{*}}C\mathrm{d}z.
```
Computing this for both salinity and temperature tracers then gives the density flux
```math
F_{\rho} = \rho_{0} \left(αF_{T} - βF_{S}\right).
```
I have taken ``α`` and ``\beta`` at the initial mean salinity and temperature as I did not save them (they can be easily computed from the simulation output if needed).
I need to check if I need to multiply by reference density and specific heat capacity in this.
Might not be needed if I frame in terms of buoyancy flux ``F_{b} = -g\left( αF_{T} - βF_{S}\right)``?

We also have have the tracer variance
```math
\mathrm{Var}(C) = \int_{V}C^{2}\mathrm{d}V
```
and tracer variance dissipation (which is the time change in the tracer variance).
"""

# ╔═╡ 7b5f9b88-5806-4e3a-92ae-3929af5ddc08
md"""
# Isothermal
"""

# ╔═╡ 33a637d5-5357-4fdf-bdb7-268cd5999f5d
begin
	iso = load("isothermal_fluxes_and_diff.jld2")
	Δz_iso = diff(iso["z"])
	replace!(iso["∂S∂z"], 0 => NaN)
	reverse!(iso["∂S∂z"], dims = 1)
	replace!(iso["κₛ"], Inf => NaN)
	replace!(iso["κₛ"], 0 => NaN)
	replace!(iso["κₛ"], -Inf => NaN)
	reverse!(iso["κₛ"], dims = 1)
	reverse!(iso["Fₜ"], dims = 1)
	reverse!(iso["Fₛ"], dims = 1)
	keys(iso)
end

# ╔═╡ 4079a130-62b6-402d-a2e7-940192874f5d
md"""
## Effecitve diffusivity

Horizontally averaging before computing the effective diffusivity leads to some *negative values* (these come from the salt flux changing sign).
In the two plots below I have taken the absolute value so we can see the magnitude to get an idea of why the effective diffusivity increases, to the parameterised molecular value, as the simulation runs.

In the volume integrated effective diffusivity there is onle **one negative value** early on in the simulation.

*There is also something happening here with the calculation of effective diffusivity.
When I do it manually and plot it I get exactly what I should but the saved computed value does not match.
Need to go through line by line and see what is going on.*
"""

# ╔═╡ 1155f907-8771-408d-b525-0d0abbb2a780
let
	fig, ax, hm = heatmap(iso["time"][1:end-1], iso["z"], log10.(iso["∂S∂z"]'), colormap = :haline)
	ax.title = "Horizontally averaged salinity gradient"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm, label = "∂S∂z (log10)")
	fig
end

# ╔═╡ a2fb3348-1272-4860-b364-7c34782ff0a5
let
	fig, ax, hm = heatmap(iso["time"][1:end-1], iso["z"], log10.(abs.(iso["κₛ"]')),
						colorrange = (log10(1e-8), log10(1)), colormap = :tempo, 
						highclip = :red, lowclip = :orange)
	ax.title = "Horizontally averaged effective diffusivity"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm, label = "Effective diffusivity (m²s⁻¹, log10)")
	fig
end

# ╔═╡ def05ee6-781c-42e4-8700-70ca78c41925
let
	# fig, ax = lines(iso["time"][1:end-1], iso["∫κₛ"])
	# hlines!(ax, 1e-7, color = :black, linestyle = :dash)
	fig, ax = lines(iso["time"][1:end-1], log10.(abs.(iso["∫κₛ"])))
	hlines!(ax, log10.(1e-7), color = :black, linestyle = :dash)
	ax.title = "Integrated horizontally averaged effective diffusivity"
	ax.xlabel = "time (s)"
	ax.ylabel = "Effective diffusivity (m²s⁻¹)"
	fig
end

# ╔═╡ fe6cade8-dd0f-4b16-aa37-babb4a2ae791
md"""
This estimate is not bad but I realised that we want the *average* effective diffusivity but the length of the domain we are averaging over **changes**.
Hence to get an accurate estimate of the average we also need to take into account the **length of the domain we are averaging over**,
```math
\kappa_{\mathrm{eff}}(t) = \frac{\int_{z \neq \texttt{NaN}}\kappa_{\mathrm{eff}}(z, t)\mathrm{d}z}{\int_{z \neq \texttt{NaN}} \mathrm{d}z}.
```

Implementing this averaging gives a much more consistent estimate of the parameterised molecular diffusivity which controls the salinity evolution in the isothermal experiment.
"""

# ╔═╡ 8cd975f4-0b08-4d5b-9acf-541e19c76391
begin
	test_κ_iso = similar(iso["∫κₛ"])
	i = 1
	for c ∈ eachcol(iso["κₛ"])
		find = findall(.!isnan.(c))
		test_κ_iso[i] = sum(c[find] .* Δz_iso[find]) / sum(Δz_iso[find])
		i += 1
	end
end

# ╔═╡ 34d8db15-2124-4c80-a3b3-aedd79e08dd0
let
	# fig, ax = lines(iso["time"][1:end-1], iso["∫κₛ"])
	# hlines!(ax, 1e-7, color = :black, linestyle = :dash)
	fig, ax = lines(iso["time"][1:end-1], log10.(abs.(test_κ_iso)))
	hlines!(ax, log10.(1e-7), color = :black, linestyle = :dash)
	ax.title = "Integrated horizontally averaged effective diffusivity"
	ax.xlabel = "time (s)"
	ax.ylabel = "Effective diffusivity (m²s⁻¹)"
	fig
end

# ╔═╡ 885c5ac5-ce91-4348-981c-b48bd39845c6
md"""
To generate a moving mean I am going to use the [`TimeSeries.jl`](https://juliastats.org/TimeSeries.jl/latest/) package.
Of course this could be written myself but I have seen this package before and have wanted to try it out.
"""

# ╔═╡ 93fb2ddb-9286-4969-845b-a27cf631a82c
begin
	being_ts = 1
	times = Time(0, being_ts, 0):Minute(1):Time(11, 0, 0)
	# data = (times = times, ∫κₛ = iso["∫κₛ"][being_ts:end])
	data = (times = times, ∫κₛ = test_κ_iso)
	iso_ts = TimeArray(data, timestamp=:times)
end

# ╔═╡ 14215a5b-ffd8-491c-848b-e51bd24645a0
begin
	window_var_iso = @bind window_iso PlutoUI.Slider(5:100, show_value=true)
	md"""
	window = $window_var_iso
	"""
end

# ╔═╡ 6b204b99-e616-42da-9c33-fd92963ff711
let
	moving_μ = moving(mean, iso_ts, window_iso)
	μ = vec(values(moving_μ))
	μ_log10 = log10.(μ)
	moving_σ = moving(std, iso_ts, window_iso)
	σ = vec(values(moving_σ))
	σ_log10 = log10.(vec(values(moving_σ)))

	fig, ax = lines(iso["time"][window_iso:end-being_ts], μ_log10)
	ax.title = "Hozizontally averaged effective diffusivity"
	ax.subtitle = "Moving mean, window = $(window_iso) mins"
	ax.xlabel = "time (s)"
	ax.ylabel = "Effective diffusivity (m²s⁻¹, log10)"
	hlines!(log10.(1e-7), color = :black, linestyle = :dash)
	# lines!(ax, iso["time"][window:end-being_ts], μ_log10 .+ σ_log10)
	# lines!(ax, iso["time"][window:end-being_ts], μ_log10 .- σ_log10)
	fig
end

# ╔═╡ 3af61b46-9a3e-49af-8d4a-ebee72f62c96
md"""
## Salinity variance dissipation
"""

# ╔═╡ 2bff8570-e034-4b2d-bd6f-0042cea0d0a0
let
	fig, ax = lines(iso["time"], log10.(iso["∫S²dV"]))
	ylims!(ax, 1.08059 .+ (-3e-6, 1e-9))
	ax.title = "Salinity variance"
	ax.xlabel = "time (s)"
	ax.ylabel = "∫S²dV"
	fig
end

# ╔═╡ 8cdf7dc6-6fe9-4c31-95ac-5841314e8c01
let
	fig, ax = lines(iso["time"][2:end], iso["dₜ∫S²dV"])
	ax.title = "Salinity variance dissipation"
	ax.xlabel = "time (s)"
	ax.ylabel = "dₜ∫S²dV"
	fig
end

# ╔═╡ c8895415-4d2a-44f7-a15f-8dd77503d940
md"""
## Density flux

This is prortional to salinity flux in this case because the temperature is isothermal over the domain.

"""

# ╔═╡ 12a475e5-4937-4e88-90b0-1d4f9a84e391
begin
	ρ₀ = 1027
	α_iso, β_iso = gsw_alpha(34.7, 0.5, 0), gsw_beta(34.7, 0.5, 0)
	F_ρ_iso = @. ρ₀ * (α_iso * iso["Fₜ"] - β_iso *  iso["Fₛ"])
	∫F_ρ_iso= vec(sum(F_ρ_iso .* Δz_iso[1], dims = 1))
	nothing
end

# ╔═╡ 5f5e360a-1dec-4f62-b45c-b771ecb0a7da
let
	F_ρ_iso_shift = log10.(F_ρ_iso .+ 1) .- log10(1)
	fig, ax, hm = heatmap(iso["time"], iso["z"], F_ρ_iso_shift')
	ax.title = "Horizontally averaged density flux"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm)
	fig
end

# ╔═╡ 483b7fb3-07e4-4c86-9c5b-15c6d437d02d
let
	fig, ax = lines(iso["time"][2:end], ∫F_ρ_iso)
	ax.title = "Volume integrated density flux"
	ax.xlabel = "time (s)"
	ax.ylabel = "∫F_ρdV"
	fig
end

# ╔═╡ 060b7142-ee7c-43cb-9892-68fd8a332daf
md"""
# Cabbeling

Need to check the salinity flux saving here --- looks like there is an error and everything is just zero which it definitely should not be.
"""

# ╔═╡ b7f3c4b1-b9b1-4dfa-8e45-5c20adf0d868
begin
	cab = load("cabbeling_fluxes_and_diff.jld2")
	Δz_cab = diff(cab["z"])
	replace!(cab["∂S∂z"], 0 => NaN)
	reverse!(cab["∂S∂z"], dims = 1)
	replace!(cab["κₛ"], Inf => NaN)
	replace!(cab["κₛ"], 0 => NaN)
	replace!(cab["κₛ"], -Inf => NaN)
	reverse!(cab["κₛ"], dims = 1)
	reverse!(cab["Fₜ"], dims = 1)
	reverse!(cab["Fₛ"], dims = 1)
	keys(cab)
end

# ╔═╡ c6dc6652-487a-4b5e-9a5f-94563e768f6e
md"""
## Effective diffusivity
"""

# ╔═╡ 885a9ae3-a19b-47c9-9eea-31d2653888ad
let
	fig, ax, hm = heatmap(cab["time"][1:end-1], cab["z"], log10.(cab["∂S∂z"]'), colormap = :haline)
	ax.title = "Horizontally averaged salinity gradient"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm, label = "∂S∂z (log10)")
	fig
end

# ╔═╡ 62353120-49c2-4a2a-9ec0-d72249fd1280
let
	fig, ax, hm = heatmap(cab["time"][1:end-1], cab["z"], log10.(abs.(cab["κₛ"]')),
						colorrange = (log10(1e-8), log10(1)), colormap = :tempo, 
						#highclip = :red, lowclip = :orange
	)
	ax.title = "Horizontally averaged effective diffusivity"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm, label = "Effective diffusivity (m²s⁻¹, log10)")
	fig
end

# ╔═╡ 75f0c7a3-04ff-4679-bd79-8ba911430c30
let
	fig, ax = lines(cab["time"][1:end-1], log10.(abs.(cab["∫κₛ"])))
	hlines!(ax, log10.(1e-7), color = :black, linestyle = :dash)
	# fig, ax = lines(cab["time"][1:end-1], cab["∫κₛ"])
	# hlines!(ax, 1e-7, color = :black, linestyle = :dash)
	ax.title = "Integrated horizontally averaged effective diffusivity"
	ax.xlabel = "time (s)"
	ax.ylabel = "Effective diffusivity (m²s⁻¹)"
	fig
end

# ╔═╡ f246bb02-7cb2-461f-9c39-9c458601ef33
begin
	test_κ_cab = similar(cab["∫κₛ"])
	j = 1
	for c ∈ eachcol(cab["κₛ"])
		find = findall(.!isnan.(c))
		test_κ_cab[j] = sum(c[find] .* Δz_cab[find]) / sum(Δz_cab[find])
		j += 1
	end
end

# ╔═╡ 3d3879a6-22ba-4107-86b3-32f6fa97050c
let
	fig, ax = lines(cab["time"][1:end-1], log10.(abs.(test_κ_cab)))
	hlines!(ax, log10.(1e-7), color = :black, linestyle = :dash)
	# fig, ax = lines(cab["time"][1:end-1], cab["∫κₛ"])
	# hlines!(ax, 1e-7, color = :black, linestyle = :dash)
	ax.title = "Integrated horizontally averaged effective diffusivity"
	ax.xlabel = "time (s)"
	ax.ylabel = "Effective diffusivity (m²s⁻¹)"
	fig
end

# ╔═╡ 757f818c-8b04-4dec-9f10-410194954384
begin
	being_ts_cab = 1
	times_cab = Time(0, being_ts_cab, 0):Minute(1):Time(11, 0, 0)
	# data_cab = (times = times_cab, ∫κₛ = cab["∫κₛ"][being_ts_cab:end])
	data_cab = (times = times_cab, ∫κₛ = test_κ_cab)
	cab_ts = TimeArray(data_cab, timestamp=:times)
end

# ╔═╡ 40f780eb-7ae5-4251-a302-a5016f13d6db
begin
	window_var_cab = @bind window_cab PlutoUI.Slider(5:100, show_value=true)
	md"""
	window = $window_var_cab
	"""
end

# ╔═╡ 3621f66f-b16b-4c6e-aa40-566f6a1fe752
let
	moving_μ = moving(mean, cab_ts, window_cab)
	μ = vec(values(moving_μ))
	μ_log10 = log10.(μ)
	moving_σ = moving(std, cab_ts, window_cab)
	σ = vec(values(moving_σ))
	σ_log10 = log10.(vec(values(moving_σ)))

	fig, ax = lines(cab["time"][window_cab:end-being_ts_cab], μ_log10)
	ax.title = "Hozizontally averaged effective diffusivity"
	ax.subtitle = "Moving mean, window = $(window_cab) mins"
	ax.xlabel = "time (s)"
	ax.ylabel = "Effective diffusivity (m²s⁻¹, log10)"
	hlines!(log10.(1e-7), color = :black, linestyle = :dash)
	# lines!(ax, cab["time"][window:end-being_ts], μ_log10 .+ σ_log10)
	# lines!(ax, cab["time"][window:end-being_ts], μ_log10 .- σ_log10)
	fig
end

# ╔═╡ a64c90c5-f75d-4ac1-84a4-cbd48fb3082b
md"""
## Tracer variance dissipation
"""

# ╔═╡ 0942924b-aa60-4946-af97-fbea014c3e9c
let
	fig, ax = lines(cab["time"], cab["∫S²dV"])
	ax.title = "Salinity variance"
	ax.xlabel = "time (s)"
	ax.ylabel = "∫S²dV"
	hidexdecorations!(ax, grid = false)
	ax2 = Axis(fig[2, 1], xlabel = "time (s)", ylabel = "∫T²dV", title = "Temperature variance")
	lines!(ax2, cab["time"], cab["∫T²dV"])
	fig
end

# ╔═╡ 77e44132-8a97-4b96-bc39-9fb8efcbd07f
let
	fig, ax = lines(cab["time"][2:end], cab["dₜ∫S²dV"])
	ax.title = "Salinity variance dissipation"
	ax.xlabel = "time (s)"
	ax.ylabel = "∫S²dV"
	hidexdecorations!(ax, grid = false)
	ax2 = Axis(fig[2, 1], xlabel = "time (s)", ylabel = "∫T²dV", title = "Temperature variance dissipation")
	lines!(ax2, cab["time"][2:end], cab["dₜ∫T²dV"])
	fig
end

# ╔═╡ 93641a6d-1aad-4651-bad3-d06297fa6342
md"""
## Density flux
"""

# ╔═╡ a74d95ff-6a44-4f60-b31e-c7a196eb7aab
begin
	g = 9.81
	Sₘ, Tₘ = 0.5 * (34.58 + 34.7), 0.5 * (-1.5 + 0.5)
	α_cab, β_cab = gsw_alpha(Sₘ, Tₘ, 0), gsw_beta(Sₘ, Tₘ, 0)
	F_ρ_cab =  @. ρ₀ * (α_cab * cab["Fₜ"] - β_cab * cab["Fₛ"])
	∫F_ρ_cab = vec(sum(F_ρ_cab .* Δz_cab[1], dims = 1))
	nothing
end

# ╔═╡ 3701b3aa-d28e-47cc-9539-d0327995250f
let
	fig, ax, hm = heatmap(iso["time"], iso["z"], F_ρ_cab')
	# fig, ax, hm = heatmap(iso["time"], iso["z"], log10.(abs.(F_ρ_cab)'))
	ax.title = "Horizontally averaged density flux"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	Colorbar(fig[1, 2], hm)
	fig
end

# ╔═╡ eef3bde2-9709-4bd3-894e-801f2db806a7
let
	fig, ax = lines(cab["time"][2:end], ∫F_ρ_cab)
	ax.title = "Volume integrated density flux"
	ax.xlabel = "time (s)"
	ax.ylabel = "∫F_ρdV"
	fig
end

# ╔═╡ 4c613a5a-cc72-4846-81f9-4fd5f9cc72e6
md"""
The volume integrated density/buoyancy flux should be related to the kinetic energy and turbulent kinetic energy via (Winters et al. (1995))
```math
\frac{\mathrm{d}}{\mathrm{d} t}E_{k} = -g\int_{V}ρw\mathrm{d}V - \epsilon.
```
The difference in magnitude means density flux is negligible compared to ``\epsilon`` so not sure if I am missing a reference density somewhere to make these things closer or it is just much smaller.
According to Inoue et al. bouyancy flux is defined as ``-g(αF_{T} - βF_{S})``.
"""

# ╔═╡ f1e195a5-ea3c-4898-9709-7bd9855bbdef
begin
	cab_energy_path = "../outputs_equaldiffusion/cabbeling_stepchange_nothing_660min/cabbeling_energetics.jld2"
	cab_energy = load(cab_energy_path)
	nothing
end

# ╔═╡ 11f01195-8bb0-4835-a722-b50e06efb314
begin

	dₜEk = diff(cab_energy["Ek"]) ./ diff(cab_energy["time"])
	ϵ = cab_energy["ϵ"]

	RHS = diff(ϵ)[1:end-1] .- (diff((g / ρ₀) * ∫F_ρ_cab) ./ diff(cab_energy["time"][1:end-1]))

	fig, ax = lines(cab_energy["time"][2:end], dₜEk, label = "dₜEk")
	lines!(ax, cab_energy["time"][2:end-1], RHS, label = "ϵ + g∫ρwdV", linestyle = :dash)
	ax.title = "Time change in kinetic energy"
	ax.xlabel = "time (s)"
	ax.ylabel = "Watts (J/s)"
	axislegend(ax)
	fig
end

# ╔═╡ cb752927-287f-4e57-b4fc-0a19777bf1e5
TableOfContents(title="Horizontally averaged fluxes and diff")

# ╔═╡ Cell order:
# ╟─b50f0e02-5ea4-11ef-0758-95c6ef87071e
# ╟─d11f7aad-5cce-4957-a9cb-8366e2b219f1
# ╟─7b5f9b88-5806-4e3a-92ae-3929af5ddc08
# ╟─33a637d5-5357-4fdf-bdb7-268cd5999f5d
# ╟─4079a130-62b6-402d-a2e7-940192874f5d
# ╟─1155f907-8771-408d-b525-0d0abbb2a780
# ╟─a2fb3348-1272-4860-b364-7c34782ff0a5
# ╟─def05ee6-781c-42e4-8700-70ca78c41925
# ╟─fe6cade8-dd0f-4b16-aa37-babb4a2ae791
# ╟─8cd975f4-0b08-4d5b-9acf-541e19c76391
# ╟─34d8db15-2124-4c80-a3b3-aedd79e08dd0
# ╟─885c5ac5-ce91-4348-981c-b48bd39845c6
# ╟─93fb2ddb-9286-4969-845b-a27cf631a82c
# ╟─14215a5b-ffd8-491c-848b-e51bd24645a0
# ╟─6b204b99-e616-42da-9c33-fd92963ff711
# ╟─3af61b46-9a3e-49af-8d4a-ebee72f62c96
# ╟─2bff8570-e034-4b2d-bd6f-0042cea0d0a0
# ╟─8cdf7dc6-6fe9-4c31-95ac-5841314e8c01
# ╟─c8895415-4d2a-44f7-a15f-8dd77503d940
# ╟─12a475e5-4937-4e88-90b0-1d4f9a84e391
# ╟─5f5e360a-1dec-4f62-b45c-b771ecb0a7da
# ╟─483b7fb3-07e4-4c86-9c5b-15c6d437d02d
# ╟─060b7142-ee7c-43cb-9892-68fd8a332daf
# ╟─b7f3c4b1-b9b1-4dfa-8e45-5c20adf0d868
# ╟─c6dc6652-487a-4b5e-9a5f-94563e768f6e
# ╟─885a9ae3-a19b-47c9-9eea-31d2653888ad
# ╠═62353120-49c2-4a2a-9ec0-d72249fd1280
# ╟─75f0c7a3-04ff-4679-bd79-8ba911430c30
# ╟─f246bb02-7cb2-461f-9c39-9c458601ef33
# ╟─3d3879a6-22ba-4107-86b3-32f6fa97050c
# ╟─757f818c-8b04-4dec-9f10-410194954384
# ╟─40f780eb-7ae5-4251-a302-a5016f13d6db
# ╟─3621f66f-b16b-4c6e-aa40-566f6a1fe752
# ╟─a64c90c5-f75d-4ac1-84a4-cbd48fb3082b
# ╟─0942924b-aa60-4946-af97-fbea014c3e9c
# ╟─77e44132-8a97-4b96-bc39-9fb8efcbd07f
# ╟─93641a6d-1aad-4651-bad3-d06297fa6342
# ╟─a74d95ff-6a44-4f60-b31e-c7a196eb7aab
# ╟─3701b3aa-d28e-47cc-9539-d0327995250f
# ╟─eef3bde2-9709-4bd3-894e-801f2db806a7
# ╟─4c613a5a-cc72-4846-81f9-4fd5f9cc72e6
# ╠═11f01195-8bb0-4835-a722-b50e06efb314
# ╟─f1e195a5-ea3c-4898-9709-7bd9855bbdef
# ╟─cb752927-287f-4e57-b4fc-0a19777bf1e5
