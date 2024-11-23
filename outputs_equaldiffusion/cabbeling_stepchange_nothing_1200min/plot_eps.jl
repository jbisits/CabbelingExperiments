using JLD2, GLMakie

file = jldopen("epsilon_maximum.jld2")
ε = file["ε_max"][2:end]
close(file)

t = eachindex(ε)

fig, ax = lines(t, log10.(ε))
ax.xlabel = "time (s)"
ax.title = "Maximum in space ε"
ax2 = Axis(fig[1, 2], xlabel = "time (s)", ylabel = "η",
           title = "Kolmogorov length")
ν = 1e-6
η = (ν^3 ./ ε).^(1/4)
lines!(ax2, t, log10.(η))
fig
