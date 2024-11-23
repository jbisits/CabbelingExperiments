using JLD2, GLMakie

file = jldopen("ε_maximum.jld2")
ε = file["ε_max"][2:end]
close(file)

t = eachindex(ε)

fig, ax = lines(t, log10.(ε))
ax.xlabel = "time (s)"
ax.title = "Maximum in space ε"
ax2 = Axis(fig[1, 2], title = "Kolmogorov lenght", 
           xlabel = "time (s)", ylabel = "η")
ν = 1e-6
η = (ν^3 ./ ε).^(1/4)
lines!(t, log10.(η))
fig
