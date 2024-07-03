# A quick look at the MCG-6-30-15 spectroscopic data

using SpectralFitting
using Plots
using Gradus
using GradusSpectralModels
using JLD2
# for JLD2 compression
using CodecZlib
# using Reflionx

include("additional-models.jl")

function read_or_make_transfer_table(path; new = false)
    if isfile(path) && !new
        @info "Reading transfer function table from file..."
        f = jldopen(path, "r")
        table = f["table"]
        close(f)
        table
    else
        @info "Creating new transfer function table..."
        d = ThinDisc(0.0, Inf)
        as = collect(range(0.0, 0.998, 10))
        # 89 fails at the moment, but don't really need full parameter space for
        # MCG-6 anyway
        thetas = collect(range(1.0, 88.0, 10))
        # this takes about 300 seconds or so; using JDL2 if it exists see
        # https://github.com/astro-group-bristol/Gradus.jl/issues/194
        table = @time make_transfer_function_table(KerrMetric, d, as, thetas)
        jldopen(path, "w"; compress = true) do f
            f["table"] = table
        end
        table
    end
end

data_path = "data"
TRANSFER_FUNCTION_PATH = joinpath(data_path, "thin-disc-transfer-table-900.jld2")

table = read_or_make_transfer_table(TRANSFER_FUNCTION_PATH);

lp_model = LineProfile(
    x -> x^-3,
    table;
    E₀ = FitParam(1.0),
    a = FitParam(0.998, lower_limit = 0.0, upper_limit = 0.998),
    θ = FitParam(30.0, lower_limit = 10.0, upper_limit = 80.0),
    rout = FitParam(400.0, frozen=true),
    # speed up model evaluation about 3 times
    quadrature_points = 13,
    n_radii = 600,
)

# comparison with Laor model
model_laor = XS_Laor(K = FitParam(1.0E-5), lineE = FitParam(6.35, frozen=true), θ = FitParam(30.0, lower_limit = 10.0, upper_limit = 80.0))

# iron line and r^-3 emissivity profile
model_gradus = AsConvolution(lp_model)(DeltaLine(K = FitParam(1.0E-5), E = FitParam(6.35, frozen=true, lower_limit=6.0, upper_limit=7.0)))

domain = collect(range(1, 10, 500))

output_laor = invokemodel(domain, model_laor)
plot(domain[1:end-1], output_laor, label="Laor")

output_gradus = 0.2 * invokemodel(domain, model_gradus)
plot!(domain[1:end-1], output_gradus, label="Gradus")
