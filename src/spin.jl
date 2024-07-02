# Fits black hole spin
# MCG--6-30-15 XMM and NuSTAR data

using SpectralFitting
using Plots
using Gradus
using GradusSpectralModels
using JLD2
# for JLD2 compression
using CodecZlib
using Reflionx

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
REFLIONX_GRID_DIR = joinpath(data_path, "reflionx/grid")

refl_table = Reflionx.parse_run(REFLIONX_GRID_DIR)
table = read_or_make_transfer_table(TRANSFER_FUNCTION_PATH);

# merged XMM data (just for the sake of the example; not for publication)
spec = joinpath(data_path, "merge_xmm_grp.pha")
back = joinpath(data_path, "merge_xmm.bak")
rmf = joinpath(data_path, "merge_xmm.rsp")
xmm = XmmData(XmmEPIC(), spec, background = back, response = rmf, ancillary = nothing)
regroup!(xmm)
drop_bad_channels!(xmm)
mask_energies!(xmm, 3.0, 10.0)
subtract_background!(xmm)
normalize!(xmm)

# merged NuSTAR data
nustar_a = NuStarData(joinpath(data_path, "merge_nustar_a_grp.pha"))
regroup!(nustar_a)
drop_bad_channels!(nustar_a)
mask_energies!(nustar_a, 3.0, 50.0)
subtract_background!(nustar_a)
normalize!(nustar_a)

nustar_b = NuStarData(joinpath(data_path, "merge_nustar_b_grp.pha"))
regroup!(nustar_b)
drop_bad_channels!(nustar_b)
mask_energies!(nustar_a, 3.0, 50.0)
subtract_background!(nustar_b)
normalize!(nustar_b)

# Fit
lp_model = LineProfile(
    x -> x^-3,
    table;
    E₀ = FitParam(1.0),
    a = FitParam(0.998, lower_limit = 0.0, upper_limit = 0.998),
    # speed up model evaluation about 3 times
    quadrature_points = 13,
    n_radii = 600,
)

# use the AutoCache wrapper to avoid re-evaluating an expensive model unnecessary
model = PowerLaw() + AsConvolution(lp_model)(ReflionxTable(K = FitParam(1e-7), refl_table))
# model = PowerLaw() + AsConvolution(lp_model)(DeltaLine(E = FitParam(6.4)))
domain = collect(range(1, 10, 500))

output = invokemodel(domain, model)
plot(domain[1:end-1], output)

# always use the ISCO
model.rin_1.frozen = true
model.rout_1.frozen = true
model.E₀_1.frozen = true

# model.logξ_1.frozen = true
# model.E_cut_1.frozen = true
# model.Γ_1.frozen = true

# make sure the datasets from the same observatory are grouped together
# else the AutoCache will trigger re-eval as the domain has changed, even through the model parameters will all be the same
datasets = FittableMultiDataset(
    xmm, nustar_a, nustar_b
)
models = FittableMultiModel((model for _ in datasets.d)...)

prob = FittingProblem(models, datasets)

bind!(prob, :a_2)
bind!(prob, :a_1)
bind!(prob, :θ_1)
bind!(prob, :Γ_1)
bind!(prob, :logξ_1)
bind!(prob, :E_cut_1)

# push!(prob.bindings, [1 => 2, 1 => 8, 2 => 2, 2 => 8, 3 => 2, 3 => 8, 4 => 2, 4 => 8, 5 => 2, 5 => 8, 6 => 2, 6 => 8, 7 => 2, 7 => 8, 8 => 2, 8 => 8, 9 => 2, 9 => 8]) # <-- this does not work (was for un-merged data sets)

result = @time fit(prob, LevenbergMarquadt(); verbose = true, x_tol = 1e-3, max_iter = 100)

begin
    p = plotresult(
        datasets.d[1],
        result[1],
        xscale = :log10,
        yscale = :log10,
        xlims = (3, 30),
        ylims = (5e-5, 2),
    )
    for i = 2:lastindex(datasets.d)
        plotresult!(datasets.d[i], result[i], residual_ylims = (-5, 5))
    end
    plot(p, legend = false)
end
