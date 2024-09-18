### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

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
	load("buoyancy_flux_interp_face.jld2")
end

# ╔═╡ ffd3a976-1234-4758-8ec2-643be78af3f7
TableOfContents(title = "Energetics analysis")

# ╔═╡ Cell order:
# ╟─943cf012-755c-11ef-3493-99307e182e96
# ╟─b2ed38aa-523a-4edf-897a-5cf8c5f0ce33
# ╠═f0ebfa39-3959-4112-8575-e6238a73bfdc
# ╟─ffd3a976-1234-4758-8ec2-643be78af3f7
