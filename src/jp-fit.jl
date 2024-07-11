using SpectralFitting
using Plots
using Gradus
using GradusSpectralModels
using JLD2
# for JLD2 compression
using CodecZlib
# using Reflionx

include("additional-models.jl")

function read_transfer_path(path)
    if isfile(path)
        @info "Reading transfer function table from file..."
        f = jldopen(path, "r")
        table = f["table"]
        close(f)

        N = 10
        a_range = collect(range(0.0, 0.998, N))
        θ_range = collect(range(20.0, 60.0, N))
        eps_range = collect(range(-0.998, 0.998, N))
        Gradus.CunninghamTransferTable((a_range, θ_range, eps_range), table)
    else
        @error "Not a file"
    end
end

data_path = "data"
table = read_transfer_path(joinpath(data_path, "johannsen_transfer_table.jld2"));

for grid in table.grids
    replace!(grid.lower_f, NaN => 0.0)
    replace!(grid.upper_f, NaN => 0.0)
end

profile = begin
    f = jldopen(joinpath(data_path, "johannsen_lamp_post.jld2"), "r")
    prof = f["lamp_post"]
    close(f)
    prof
end;


model = JohannsenPsaltisLampPost(profile, table)

model.a.value = 0.6
model.eps.value = 0.0

domain = collect(range(0.1, 1.6, 200))
output = invokemodel(domain, model)

plot!(domain[1:end-1], output)
