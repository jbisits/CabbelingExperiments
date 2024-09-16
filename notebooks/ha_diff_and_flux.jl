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
	using JLD2, CairoMakie, PlutoUI, StatsBase, TimeSeries, Dates, GibbsSeaWater, Statistics
end

# ╔═╡ d11f7aad-5cce-4957-a9cb-8366e2b219f1
md"""
# Horizontally averaged effective diffusivity and fluxes

One way to get rid of the spike in the effective diffusivity is to first horizontally average.
This gives behaviour like we would expect and perhaps taking a moving mean we can generate a good estimate using this diagnostic.

The effective diffusivity at each depth level of the *horizontally averaged salinity profile* is
```math
\kappa_{\mathrm{eff}}(z^{*}, t) = \frac{∂_{t}\int_{z < z^{*}} \langle S \rangle_{xy}\mathrm{d}z}{\partial \langle S \rangle_{xy} / \partial z}.
```
Integrating over all depths and dividing by total depth, which is 1m, to get a an average effective diffusivity for a *horizontally averaged profile*.
**Note:** the horizontally averaged profile is still **sorted**.

Tracer fluxes from horizontally averaged, sorted profiles are
```math
F_{C} = ∂_{t}\int_{z < z^{*}}\langle C \rangle_{xy}\mathrm{d}z.
```
Computing this for both salinity and temperature tracers then gives the density flux
```math
F_{\rho} = \rho_{0} \left(-αF_{T} + βF_{S}\right).
```
In this calculation, the salinity/temperature field are multiplied by ``β/α`` which are computed from the *horizontally averaged tracer profiles*.
Then, to get an integrated quantity, we need to integrate over `z` and multiply by surface area,
```math
\mathrm{SA} ∫_{z}F_{\rho}\mathrm{d}z.
```
This ensures our dimensions match (I think).
I have been doing some sanity checks to see if things line up as they should.

We also have have the tracer variance
```math
\mathrm{Var}(C) = \int_{V}(C - \langle C \rangle)^{2}\mathrm{d}V
```
where ``\langle C \rangle`` is the mean of tracer ``C`` and tracer variance dissipation (which is the time change in the tracer variance).
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

# ╔═╡ 1d0dcc31-ad71-4522-8165-73194847d0c6
md"""
## Tracer fluxes
"""

# ╔═╡ f12e37c2-b18d-4684-81ec-97b10eb17138
let
	fig, ax, hm = heatmap(iso["time"][2:end], iso["z"], iso["Fₜ"]', colormap = :thermal)
	ax.title = "Temperature flux"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	hidexdecorations!(ax, ticks = false)
	Colorbar(fig[1, 2], hm)
	ax2 = Axis(fig[2, 1])
	hm = heatmap!(ax2, iso["time"][2:end], iso["z"], log10.(abs.(iso["Fₛ"]')), colormap = :haline)
	ax2.title = "Salinity flux"
	ax2.xlabel = "time (s)"
	ax2.ylabel = "z (m)"
	Colorbar(fig[2, 2], hm)
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

# ╔═╡ c6dc6652-487a-4b5e-9a5f-94563e768f6e
md"""
## Effective diffusivity
"""

# ╔═╡ 40f780eb-7ae5-4251-a302-a5016f13d6db
begin
	window_var_cab = @bind window_cab PlutoUI.Slider(5:100, show_value=true)
	md"""
	window = $window_var_cab
	"""
end

# ╔═╡ a64c90c5-f75d-4ac1-84a4-cbd48fb3082b
md"""
## Tracer variance dissipation
"""

# ╔═╡ 93641a6d-1aad-4651-bad3-d06297fa6342
md"""
## Tracer fluxes
"""

# ╔═╡ 45d9e940-39eb-408c-ba16-e923f52d1874
md"""
## Density flux
"""

# ╔═╡ 4c613a5a-cc72-4846-81f9-4fd5f9cc72e6
md"""
The volume integrated density/buoyancy flux should be related to the kinetic energy and turbulent kinetic energy via (Winters et al. (1995))
```math
\frac{\mathrm{d}}{\mathrm{d} t}E_{k} = -g\int_{V}ρw\mathrm{d}V - \epsilon.
```
"""

# ╔═╡ 0a184e81-2d7a-483c-b223-adbac5aaa234
md"""

# Sanity checks

These are dons for the cabbeling experiment.
I was having trouble figuring out closing the energy budget as per
```math
\frac{\mathrm{d}}{\mathrm{d} t}E_{k} = -g\int_{V}ρw\mathrm{d}V - \epsilon.
```
so I have computed everything again (from the in situ model output) to check and all looks to check out.
"""

# ╔═╡ d252f33a-37d8-4f68-a4bd-0a4cb58b214e
begin
	cab_long_energetics = load("buoyancy_flux.jld2")
	∫Eₖ = cab_long_energetics["∫Eₖ"]
	∫ϵ = cab_long_energetics["∫ϵ"]
	∫gρw = cab_long_energetics["∫gρw"]
	∫αΘw = cab_long_energetics["∫αΘw"]
	∫βSw = cab_long_energetics["∫βSw"]
	# ∫gρw = load("buoyancy_flux_interp_face.jld2", "∫gρw")
	ρ₀_string = cab_long_energetics["ρ₀"]
	find_num = findfirst('k', cab_long_energetics["ρ₀"]) - 1
	ρ₀_model = parse(Float64, cab_long_energetics["ρ₀"][1:find_num])
	keys(cab_long_energetics)
end

# ╔═╡ b7f3c4b1-b9b1-4dfa-8e45-5c20adf0d868
begin
	cab = load("cabbeling_fluxes_and_diff_longer_run.jld2")
	# cab = load("cabbeling_fluxes_and_diff.jld2")
	cab_time = cab["time"]
	Δz_cab = diff(cab["z"])
	replace!(cab["∂S∂z"], 0 => NaN)
	reverse!(cab["∂S∂z"], dims = 1)
	replace!(cab["κₛ"], Inf => NaN)
	replace!(cab["κₛ"], 0 => NaN)
	replace!(cab["κₛ"], -Inf => NaN)
	reverse!(cab["κₛ"], dims = 1)
	reverse!(cab["Fₜ"], dims = 1)
	reverse!(cab["Fₛ"], dims = 1)
	Cₚ = GibbsSeaWater.gsw_cp0
	interp_∫gρw = (cab_long_energetics["∫gρw"][1:end-1] .+ cab_long_energetics["∫gρw"][2:end]) / 2
	keys(cab)
end

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

# ╔═╡ 757f818c-8b04-4dec-9f10-410194954384
begin
	being_ts_cab = 1
	# times_cab = Time(0, being_ts_cab, 0):Minute(1):Time(11, 0, 0)
	times_cab = Time(0, being_ts_cab, 0):Minute(1):Time(20, 0, 0)
	# data_cab = (times = times_cab, ∫κₛ = cab["∫κₛ"][being_ts_cab:end])
	data_cab = (times = times_cab, ∫κₛ = test_κ_cab)
	cab_ts = TimeArray(data_cab, timestamp=:times)
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

# ╔═╡ 0942924b-aa60-4946-af97-fbea014c3e9c
let
	fig, ax = lines(cab["time"], cab["∫(S - Smean)²dV"])
	ax.title = "Salinity variance"
	ax.xlabel = "time (s)"
	ax.ylabel = "∫S²dV"
	hidexdecorations!(ax, grid = false)
	ax2 = Axis(fig[2, 1], xlabel = "time (s)", ylabel = "∫T²dV", title = "Temperature variance")
	lines!(ax2, cab["time"], cab["∫(T - Tmean)²dV"])
	fig
end

# ╔═╡ 77e44132-8a97-4b96-bc39-9fb8efcbd07f
let
	S_var = diff(cab["∫(S - Smean)²dV"]) ./ diff(cab_time)
	T_var = diff(cab["∫(T - Tmean)²dV"]) ./ diff(cab_time)
	fig, ax = lines(cab["time"][2:end], S_var)
	ax.title = "Salinity variance dissipation"
	ax.xlabel = "time (s)"
	ax.ylabel = "∫S²dV"
	hidexdecorations!(ax, grid = false)
	ax2 = Axis(fig[2, 1], xlabel = "time (s)", ylabel = "∫T²dV", title = "Temperature variance dissipation")
	lines!(ax2, cab["time"][2:end], T_var)
	fig
end

# ╔═╡ a78a42a5-66d1-4aad-982d-32fccb1ebb84
let
	fig, ax, hm = heatmap(cab["time"][2:end], cab["z"], cab["Fₜ"]', colormap = :thermal)
	ax.title = "Temperature flux"
	ax.xlabel = "time (s)"
	ax.ylabel = "z (m)"
	hidexdecorations!(ax, ticks = false)
	Colorbar(fig[1, 2], hm)
	ax2 = Axis(fig[2, 1])
	hm = heatmap!(ax2, cab["time"][2:end], cab["z"], cab["Fₛ"]', colormap = :haline)
	ax2.title = "Salinity flux"
	ax2.xlabel = "time (s)"
	ax2.ylabel = "z (m)"
	Colorbar(fig[2, 2], hm)
	fig
end

# ╔═╡ a74d95ff-6a44-4f60-b31e-c7a196eb7aab
begin
	Sₘ, Tₘ = 0.5 * (34.58 + 34.7), 0.5 * (-1.5 + 0.5)
	α_cab, β_cab = gsw_alpha(Sₘ, Tₘ, 0), gsw_beta(Sₘ, Tₘ, 0)
	F_ρ_cab =  ρ₀_model * (-cab["α"] .* cab["Fₜ"] .+ (cab["β"] .* cab["Fₛ"]))
	∫F_ρ_cab = vec(sum(F_ρ_cab, dims = 1)) .* Δz_cab[1] .* (0.1^2) # SA = 0.1^2
	nothing
end

# ╔═╡ 3701b3aa-d28e-47cc-9539-d0327995250f
let
	fig, ax, hm = heatmap(iso["time"], iso["z"], F_ρ_cab', colormap = :balance)
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

# ╔═╡ 11f01195-8bb0-4835-a722-b50e06efb314
begin
	
	Ek = cab_long_energetics["∫Eₖ"]
	dₜEk = diff(Ek ) ./ diff(cab_long_energetics["time"])
	ϵ =  cab_long_energetics["∫ϵ"]
	g = 9.81
	∫J_b = g * ∫F_ρ_cab / ρ₀_model

	RHS = -0.5 * (ϵ[1:end-1] .+ ϵ[2:end]) .- ∫J_b

	fig, ax = lines(cab_long_energetics["time"][1:200], dₜEk[1:200], label = "dₜEk")
	lines!(ax, cab_long_energetics["time"][1:200], RHS[1:200], label = "ϵ - ∫J_b", linestyle = :dash)
	ax.title = "Time change in kinetic energy"
	ax.xlabel = "time (s)"
	ax.ylabel = "Watts (J/s)"
	axislegend(ax)
	fig
end

# ╔═╡ 20c14d5c-b2d7-4238-a004-272eeb72014c
let
	time = cab_long_energetics["time"]
	Δt = diff(time)
	∫Eₖ = cab_long_energetics["∫Eₖ"]
	∫ϵ = cab_long_energetics["∫ϵ"]
	∫gρw = cab_long_energetics["∫gρw"]

	dₜ∫Eₖ = diff(∫Eₖ) ./ Δt
	RHS_ = -∫ϵ .- ∫gρw / ρ₀_model
	RHS = (RHS_[1:end-1] .+ RHS_[2:end]) / 2

	fig = Figure(size = (1000, 800))
	ax = Axis(fig[1, 1], title = "Kinetic energy, TKE dissipation and buoyancy flux", xlabel = "time (s)", ylabel = "Watts (J/s)")
	lines!(ax, time[1:199], dₜ∫Eₖ[1:199], label = "dₜ∫Eₖ")
	lines!(ax, time[1:199], RHS[1:199], label = "-∫ϵ - ∫gρw")
	hidexdecorations!(ax, ticks = false, grid = false)
	axislegend(ax)

	RMS = sqrt(mean((dₜ∫Eₖ[1:199] .- RHS[1:199]).^2))
	abs_err = abs.(dₜ∫Eₖ[1:199] .- RHS[1:199])
	rel_err = abs_err ./ dₜ∫Eₖ[1:199]
	ax2 = Axis(fig[2, 1], title = "Absolute error", subtitle = "RMS = $(RMS)", xlabel = "time (s)", ylabel = "Absolute error")
	lines!(ax2, time[1:199], abs_err)
	fig
end

# ╔═╡ 46fdeae9-293e-46bf-9ebc-243da60df30a
md"""
Computing from model output, so I compute ``\int_{V}ρ w \mathrm{d}V`` from the `ρ` and `w` saved output, we have good agreement between the volume integrated quantities above.

After scaling by the surface area actually have reasonably good agreement between in situ density flux and computed from fluxes density flux.
"""

# ╔═╡ 3465bbc1-b6d1-4aa9-a87f-4ed204dc9adb
let
	fig, ax = lines(cab_time[1:200], interp_∫gρw[1:200] / cab_long_energetics["g"], label = "∫ρw")
	lines!(ax, cab_time[1:200], ∫F_ρ_cab[1:200], label = "From sorted S, T fluxes")
	# lines!(ax, cab_time[1:200], -F_ρ_vol_integrated[1:200])
	ax.title = "Comparison between density flux from S and T fluxes and computed directly from model"
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ 5d65a2c1-c535-43e0-92aa-e447c979ae9c
md"""
If we compute ``\int_{V}α Θ w\mathrm{d}V``, ``\int_{V}β S w\mathrm{d}V``, from the in-situ model output then
```math
ρ₀\left(-∫_{V}αΘw\mathrm{d}V + ∫_{V} βSw\mathrm{d}V\right)
```
we should also have good agreement, but do not.

There are still some dimension things going on that I have not sorted.
"""

# ╔═╡ 032c0b8b-1363-4966-a963-80c141f5de6c
let
	F_ρ_vol_integrated = -∫αΘw .+ ∫βSw
	interp_F_ρ_vol_integrated = ρ₀_model * 0.5 * (F_ρ_vol_integrated[1:end-1] .+ F_ρ_vol_integrated[2:end])
	
	fig, ax = lines(cab_long_energetics["time"][1:200], ∫F_ρ_cab[1:200], label = "From sorted S, T fluxes")
	lines!(ax, cab_long_energetics["time"][1:200], interp_F_ρ_vol_integrated[1:200], label = "ρ₀(∫αΘw - ∫βSw)")
	ax.xlabel = "time (s)"
	ax.ylabel = ""
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ 6162b2c8-bcba-466d-90fe-4911a6e17122
md"""
Not seeing these lineup even though now the density flux looks to be reasnable!

Try computing only using salinity so I compare ``\mathrm{SA} \int_{z}ρ_{0}\left(βF_{S} \right) \mathrm{d}z`` to ``ρ₀\int_{V}βSw\mathrm{d}V``.
The signal here is matching well but underestimating using the horizontally averaged salinity profile and flux computed from that.

The difference is coming from the salinity more than the temperature but this difference does not account for the difference from the total density flux in the above figure.

After computing with these quantities seperately all seems to work out after all so there must be an error above.

This all boils down to **where negative signs are**.
I need to match this up and get a dimensional diagnosis from Jan on this.

The integrated temperature fluxes have the *same sign* but the salt fluxes have *alternate signs*.
This could be where the issue is.
The odd thing is the volume integrated tracers are out by a factor of about 15...
"""

# ╔═╡ 867fe3cf-c1d6-4470-bfd6-f3381de73f21
let
	ρ_S_fluxes = 0.1^2 * ρ₀_model * sum(cab["β"] .* cab["Fₛ"], dims = 1) * Δz_cab[1]
	ρ_S_vol_int = -0.5 * (ρ₀_model * (∫βSw[1:end-1] + ∫βSw[2:end]))
	
	fig, ax = lines(cab["time"][1:200], vec(ρ_S_fluxes)[1:200], label = "SAρ₀∫βFₛdz")
	lines!(ax, cab["time"][1:200], ρ_S_vol_int[1:200], label = "-ρ₀∫βSwdV")
	ax.title = "''Density'' from salinity"
	ax.xlabel = "time (s)"
	ax.ylabel = ""
	axislegend(ax, position = :rt)
	
	ρ_T_fluxes = -0.1^2 * ρ₀_model * sum(cab["α"] .* cab["Fₜ"], dims = 1) * Δz_cab[1]
	ρ_T_vol_int = -0.5 * (ρ₀_model * (∫αΘw[1:end-1] + ∫αΘw[2:end]))

	ax2 = Axis(fig[2, 1], xlabel = "time (s)", title = "''Density'' from temperature")
	lines!(ax2, cab["time"][1:200], vec(ρ_T_fluxes)[1:200], label = "-SAρ₀∫αFₜdz")
	lines!(ax2, cab["time"][1:200], ρ_T_vol_int[1:200], label = "-ρ₀∫αΘwdV")
	axislegend(ax2, position = :rb)

	ρ_fluxes = ρ_T_fluxes .+ ρ_S_fluxes
	ρvol_int = (ρ_T_vol_int .- ρ_S_vol_int)

	ax3 = Axis(fig[:, 2], xlabel = "time (s)", title = "Density")
	lines!(ax3, cab["time"][1:200], vec(ρ_fluxes)[1:200], label = "SAρ₀∫(αFₜ - βFₛ)dz")
	lines!(ax3, cab["time"][1:200], interp_∫gρw[1:200] / cab_long_energetics["g"] ,label = "∫ρwdV")
	lines!(ax3, cab["time"][1:200], ρvol_int[1:200], label = "ρ₀∫αΘwdV - ρ₀∫βSwdV")
	axislegend(ax3, position = :rb)
	fig
end

# ╔═╡ 4b5d1c49-ca72-4f1a-9737-7a1f8d4411c9
md"""
## Energy budget

I now have everything to close the energy budget.
First I check a few different ways to compute the energetics to make sure things are matching up.

### Potential energy

One version is computed from a `KernelFunctionOperation` I wrote for [`Oceanostics.jl`](https://github.com/tomchor/Oceanostics.jl/blob/51f477ed1e3d9d2ab0683c47609f2f8e2afd6243/src/PotentialEnergyEquationTerms.jl#L180) where I compute things on the fly and another is based off the saved output where I compute
```math
	E_{p} = g∫_{V}ρz\mathrm{d}V \tag{1}
```
I use potential density ``σ_{0}``.
`∫Ep_alternate` is the method computed from the saved output using ``(1)``.
From figures below can see there is a fairly consistent difference between the two --- which one is correct?

More specifically which one should I use?
Perhaps it does not matter as the time changing quantities are almost exactly equal.
"""

# ╔═╡ 61770c4f-f8ca-46da-9ff0-73d5aadd2972
begin
	∫Ep = load("cabbeling_fluxes_and_diff_longer_run.jld2", "∫Ep")
	∫Ep_alternate = load("cabbeling_fluxes_and_diff_longer_run.jld2", "∫Ep_alternate")
	t = load("cabbeling_fluxes_and_diff_longer_run.jld2", "time")
	t_interp = 0.5 * (t[1:end-1] .+ t[2:end])
	dₜ∫Ep = diff(∫Ep) ./ diff(t)
	dₜ∫Ep_alternate = diff(∫Ep_alternate) ./ diff(t)
	dₜ∫Ek = diff(∫Eₖ) ./ diff(t)
	nothing
end

# ╔═╡ c852c48d-a715-4a37-8456-e5a18fccf47b
let
	fig, ax = lines(t, ∫Ep, label = "Computed on fly")
	ax.ylabel = "Joules"
	ax.title = "Potential energy calculations"
	hidexdecorations!(ax, grid=false, ticks=false)
	axislegend(ax)
	ax2 = Axis(fig[2, 1], ylabel = "Joules")
	lines!(ax2, t, ∫Ep_alternate, label = "Computed post", color = :orange)
	hidexdecorations!(ax2, grid=false, ticks=false)
	axislegend(ax2)
	RMSE = sqrt(mean(∫Ep .- ∫Ep_alternate).^2)
	abs_err = abs.(∫Ep .- ∫Ep_alternate)
	ax.subtitle = "RMSE = $(RMSE)"
	ax3 = Axis(fig[3, 1], ylabel = "absolute error")
	lines!(ax3, t, abs_err)
	ax3.xlabel = "time (s)"
	fig
end

# ╔═╡ 101a2635-1198-46b4-9595-87de587919fa
let
	fig, ax = lines(t_interp, dₜ∫Ep, label = "Computed on fly")
	ax.ylabel = "Watts"
	ax.title = "Time change potential energy calculations"
	hidexdecorations!(ax, grid=false, ticks=false)
	lines!(ax, t_interp, dₜ∫Ep_alternate, label = "Computed post", linestyle = :dash)
	axislegend(ax)
	RMSE = sqrt(mean(dₜ∫Ep .- dₜ∫Ep_alternate).^2)
	abs_err = abs.(dₜ∫Ep .- dₜ∫Ep_alternate)
	ax.subtitle = "RMSE = $(RMSE)"
	ax2 = Axis(fig[2, 1], ylabel = "absolute error")
	lines!(ax2, t_interp, abs_err)
	ax2.xlabel = "time (s)"
	fig
end

# ╔═╡ 6f3dc607-114a-4f7a-a51d-643fa8ad8fa6
md"""
### Background potential energy

Again I compute this two ways.
One way is defining
```math
z✶ = \frac{1}{SA}∫_{ρ_{\mathrm{min}}}^{ρ_{\mathrm{max}}}\mathrm{d}V
```
then
```math
E_{b} = g∫_{V}ρ✶z✶\mathrm{d}V
```
where ``ρ✶`` is the reshaped and sorted density profile.
The other method I use is to create a vector of heights so that each element in the reshaped one d density profile is associated to a height.
Then after the profile is sorted integrate over all depth.

They end up with opposite signs but mathcing them have good agreement.

**Currently the sign switching is actually a problem** I need to make usre I compute the BPE and the PE in the same manner ---  that is, in the APE calculation I have ``z✶`` going from 0...1, so densest water is at 0 and lightest water is at one, but in the potential energy case I compute with densest water at -1.
This is why what I had as the previous version of BPE (with `rev = false` in `sort!`) worked but this current one, where I compute what I actually think I should is not correct.
I think this boils down to offsetting the grid on which the PE is calculated.
This will be better as then will see APE rising as time goes on (which is what the plot shows).
"""

# ╔═╡ 85773eb5-d727-4a73-88e7-ec0ee6c7eda2
begin
	∫Ebz✶ = -load("cabbeling_fluxes_and_diff_longer_run.jld2", "∫Ebz✶")
	∫Eb = load("cabbeling_fluxes_and_diff_longer_run.jld2", "∫Eb") ./ ρ₀_model
	dₜ∫Ebz✶ = diff(∫Ebz✶) ./ diff(t)
	dₜ∫Eb = diff(∫Eb) ./ diff(t)
	nothing
end

# ╔═╡ 1293f605-90ab-4f18-934c-4718a8eec6a9
let
	fig, ax = lines(t, ∫Eb, label = "∫Eb")
	ax.ylabel = "Joules"
	ax.title = "Background potential energy calculations"
	hidexdecorations!(ax, grid=false, ticks=false)
	lines!(ax, t, ∫Ebz✶, label = "∫Ebz✶", color = :orange)
	axislegend(ax)
	RMSE = sqrt(mean(∫Eb .- ∫Ebz✶).^2)
	abs_err = abs.(∫Eb .- ∫Ebz✶)
	ax.subtitle = "RMSE = $(RMSE)"
	ax2 = Axis(fig[2, 1], ylabel = "absolute error")
	lines!(ax2, t, abs_err)
	ax2.xlabel = "time (s)"
	fig
end

# ╔═╡ 467ed030-3e4f-4d80-a593-5205706ce633
let
	fig, ax = lines(t_interp[1:200], dₜ∫Eb[1:200], label = "dₜ∫Eb")
	ax.ylabel = "Watts"
	ax.title = "Background potential energy calculations"
	hidexdecorations!(ax, grid=false, ticks=false)
	lines!(ax, t_interp[1:200], dₜ∫Ebz✶[1:200], label = "dₜ∫Ebz✶", linestyle = :dash)
	axislegend(ax)
	RMSE = sqrt(mean(dₜ∫Eb .- dₜ∫Ebz✶).^2)
	abs_err = abs.(dₜ∫Eb .- dₜ∫Ebz✶)
	ax.subtitle = "RMSE = $(RMSE)"
	ax2 = Axis(fig[2, 1], ylabel = "absolute error")
	lines!(ax2, t_interp, abs_err)
	ax2.xlabel = "time (s)"
	fig
end

# ╔═╡ 7d9c1960-844d-4362-aa56-d4790fb5d908
begin
	energy_window = @bind plot_window PlutoUI.Slider(2:length(∫Ep)-1, default=200)
	pe_method = @bind pe_comp PlutoUI.Select(["On fly", "Post processing"])
	bpe_method = @bind bpe_comp PlutoUI.Select(["z✶", "repeated_z"])
	md"""
	### Available potential energy
	
	Window over which to view energetics.
	$(energy_window)
	
	Method used for ``\frac{\mathrm{d}}{\mathrm{d}t}E_{p}``: $(pe_method).
	
	Method used for ``\frac{\mathrm{d}}{\mathrm{d}t}E_{b}``: $(bpe_method).
	"""
end

# ╔═╡ 05fb4be9-c1f8-415c-b5a8-ceaf8853d818
let
	plot_dₜEb = bpe_comp == "z✶" ? dₜ∫Ebz✶ : dₜ∫Eb
	plot_dₜEp = pe_comp == "On fly" ? dₜ∫Ep : dₜ∫Ep_alternate
	plot_dₜEa = plot_dₜEp .- plot_dₜEb
	fig = Figure(size = (500, 900)) 
	ax = Axis(fig[1, 1], xlabel = "time (s)", ylabel = "Watts", title = "Time change in potential energies", subtitle = "window = 0-$(t[plot_window] / 60)min")
	lines!(ax, t_interp[1:plot_window], plot_dₜEp[1:plot_window], label = "dₜ∫Ep")
	lines!(ax, t_interp[1:plot_window], plot_dₜEb[1:plot_window], label = "dₜ∫Eb")
	hidexdecorations!(ax, grid = false, ticks = false)
	axislegend(ax, position = :rb)
	ax2 = Axis(fig[2, 1], xlabel = "time (s)", ylabel = "Watts", title = "Time change in APE", subtitle = "window = 0-$(t[plot_window] / 60)min")
	lines!(ax2, t_interp[1:plot_window], plot_dₜEa[1:plot_window], color = :red, label = "dₜ∫Ea")
	lines!(ax2, t_interp[1:plot_window], dₜ∫Ek[1:plot_window], color = :green, label = "dₜ∫Ek")
	axislegend(ax2)
	fig
end

# ╔═╡ 5f322ee6-f974-4005-9acf-7a8078b2eca3
begin
	Φd = dₜEb = bpe_comp == "z✶" ? dₜ∫Ebz✶ : dₜ∫Eb
	dₜEp = pe_comp == "On fly" ? dₜ∫Ep : dₜ∫Ep_alternate
	dₜEa = dₜEp .- dₜEb
	∫ϵ_interp = 0.5 * (∫ϵ[1:end-1] .+ ∫ϵ[2:end])
	Φz = interp_∫gρw ./ ρ₀
	Φi = dₜ∫Ep .- Φz 
end

# ╔═╡ 4a328b42-45f6-4a03-87dc-6e8efd70a83d
let
	fig = Figure(size = (500, 500))
	ax = Axis(fig[1, 1], xlabel = "time (s)", ylabel = "Watts", title = "Fluxes")
	lines!(ax, t_interp[1:200], Φz[1:200], label = "Φz")
	lines!(ax, t_interp[1:200], Φd[1:200], label = "Φd")
	lines!(ax, t_interp[1:200], ∫ϵ_interp[1:200], label = "∫ϵ")
	lines!(ax, t_interp[1:200], Φi[1:200], label = "Φi")
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ d97f661d-7133-423c-9dda-81bad26659e0
md"""
## Background potential energy from salinity and temperature

The change in the background potential energy is due to a change in distribution of the density.
In Winters (1995) they have linear (i.e. conserved) density so the only way this distribution changes is due to irreversible mixing.

I want to know how the density distribution changes due linear and non-linear processes.
A way to estimate the non-linear effect is to *subtract the linear change in BPE from the total BPE*,
```math
\begin{aligned}
\frac{\mathrm{d}}{\mathrm{d} t}E_{b} &= \frac{\mathrm{d}}{\mathrm{d} t}E_{b}^{\mathrm{linear}} + \frac{\mathrm{d}}{\mathrm{d} t}E_{b}^{\mathrm{non-linear}} \\
\implies \frac{\mathrm{d}}{\mathrm{d} t}E_{b}^{\mathrm{non-linear}} &= \frac{\mathrm{d}}{\mathrm{d} t}E_{b} -\frac{\mathrm{d}}{\mathrm{d} t}E_{b}^{\mathrm{linear}}
\end{aligned}
```

I have the total BPE so just need the linear change in BPE (not sure if this is good terminology).
"""

# ╔═╡ 8d5c31ab-c16c-433e-bf74-655423595b84
let
	S_flux = 0.5 * (∫βSw[1:end-1] + ∫βSw[2:end])
	T_flux = -0.5 * (∫αΘw[1:end-1] + ∫αΘw[2:end])
	S_and_T_fluxes = S_flux .+ T_flux
	fig, ax = lines(t_interp[1:200], S_and_T_fluxes[1:200], label = "S and T fluxes")
	lines!(ax, t_interp[1:200], dₜ∫Ebz✶[1:200], label = "BPE")
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ f1e195a5-ea3c-4898-9709-7bd9855bbdef
begin
	cab_energy_path = "../outputs_equaldiffusion/cabbeling_stepchange_nothing_660min/cabbeling_energetics.jld2"
	cab_energy = load(cab_energy_path)
	nothing
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
# ╟─1d0dcc31-ad71-4522-8165-73194847d0c6
# ╟─f12e37c2-b18d-4684-81ec-97b10eb17138
# ╟─c8895415-4d2a-44f7-a15f-8dd77503d940
# ╟─12a475e5-4937-4e88-90b0-1d4f9a84e391
# ╟─5f5e360a-1dec-4f62-b45c-b771ecb0a7da
# ╟─483b7fb3-07e4-4c86-9c5b-15c6d437d02d
# ╟─060b7142-ee7c-43cb-9892-68fd8a332daf
# ╟─b7f3c4b1-b9b1-4dfa-8e45-5c20adf0d868
# ╟─c6dc6652-487a-4b5e-9a5f-94563e768f6e
# ╟─885a9ae3-a19b-47c9-9eea-31d2653888ad
# ╟─62353120-49c2-4a2a-9ec0-d72249fd1280
# ╟─f246bb02-7cb2-461f-9c39-9c458601ef33
# ╟─3d3879a6-22ba-4107-86b3-32f6fa97050c
# ╟─757f818c-8b04-4dec-9f10-410194954384
# ╟─40f780eb-7ae5-4251-a302-a5016f13d6db
# ╟─3621f66f-b16b-4c6e-aa40-566f6a1fe752
# ╟─a64c90c5-f75d-4ac1-84a4-cbd48fb3082b
# ╟─0942924b-aa60-4946-af97-fbea014c3e9c
# ╟─77e44132-8a97-4b96-bc39-9fb8efcbd07f
# ╟─93641a6d-1aad-4651-bad3-d06297fa6342
# ╟─a78a42a5-66d1-4aad-982d-32fccb1ebb84
# ╟─45d9e940-39eb-408c-ba16-e923f52d1874
# ╟─a74d95ff-6a44-4f60-b31e-c7a196eb7aab
# ╟─3701b3aa-d28e-47cc-9539-d0327995250f
# ╟─eef3bde2-9709-4bd3-894e-801f2db806a7
# ╟─4c613a5a-cc72-4846-81f9-4fd5f9cc72e6
# ╟─11f01195-8bb0-4835-a722-b50e06efb314
# ╟─0a184e81-2d7a-483c-b223-adbac5aaa234
# ╟─d252f33a-37d8-4f68-a4bd-0a4cb58b214e
# ╟─20c14d5c-b2d7-4238-a004-272eeb72014c
# ╟─46fdeae9-293e-46bf-9ebc-243da60df30a
# ╟─3465bbc1-b6d1-4aa9-a87f-4ed204dc9adb
# ╟─5d65a2c1-c535-43e0-92aa-e447c979ae9c
# ╟─032c0b8b-1363-4966-a963-80c141f5de6c
# ╟─6162b2c8-bcba-466d-90fe-4911a6e17122
# ╟─867fe3cf-c1d6-4470-bfd6-f3381de73f21
# ╟─4b5d1c49-ca72-4f1a-9737-7a1f8d4411c9
# ╟─61770c4f-f8ca-46da-9ff0-73d5aadd2972
# ╟─c852c48d-a715-4a37-8456-e5a18fccf47b
# ╟─101a2635-1198-46b4-9595-87de587919fa
# ╟─6f3dc607-114a-4f7a-a51d-643fa8ad8fa6
# ╟─85773eb5-d727-4a73-88e7-ec0ee6c7eda2
# ╟─1293f605-90ab-4f18-934c-4718a8eec6a9
# ╟─467ed030-3e4f-4d80-a593-5205706ce633
# ╟─7d9c1960-844d-4362-aa56-d4790fb5d908
# ╟─05fb4be9-c1f8-415c-b5a8-ceaf8853d818
# ╟─5f322ee6-f974-4005-9acf-7a8078b2eca3
# ╟─4a328b42-45f6-4a03-87dc-6e8efd70a83d
# ╟─d97f661d-7133-423c-9dda-81bad26659e0
# ╟─8d5c31ab-c16c-433e-bf74-655423595b84
# ╟─f1e195a5-ea3c-4898-9709-7bd9855bbdef
# ╟─cb752927-287f-4e57-b4fc-0a19777bf1e5
