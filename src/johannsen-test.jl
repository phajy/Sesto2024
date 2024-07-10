using Gradus, Plots
using SpectralFitting
include("additional-models.jl")


model = JPLineProfile(a = FitParam(0.4), θ = FitParam(30.0), a13 = FitParam(0.0), eps3 = FitParam(2.0))
domain = collect(range(0.1, 1.6, 200))
output = invokemodel(domain, model)
plot(domain[1:end-1], output)

model.eps3.value = -2

output = invokemodel(domain, model)
plot!(domain[1:end-1], output)
