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
using CSV
using DataFrames

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

model = PowerLaw()
# model = XS_CutOffPowerLaw(Ecut = FitParam(125.0, lower_limit=100.0, upper_limit=200.0), z = FitParam(0.007749, frozen=true))
# use the AutoCache wrapper to avoid re-evaluating an expensive model unnecessary
# model = PowerLaw() + AsConvolution(lp_model)(ReflionxTable(K = FitParam(1e-7), refl_table))
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

# cut off power law
# bind!(prob, :Γ)
# bind!(prob, :Ecut)
# power law
bind!(prob, :a)

# reflection model
# bind!(prob, :a_2)
# bind!(prob, :a_1)
# bind!(prob, :θ_1)
# bind!(prob, :Γ_1)
# bind!(prob, :logξ_1)
# bind!(prob, :E_cut_1)

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

# create data files for plots that can be rendered in Veusz for presentation
# note that the "pm" column names should be renamed "+-" for Veusz to interpret them as error bars
# this is a bit long winded and should be made more general
function output_result(ds, result, bin_factor, file_path)
    ds_x = SpectralFitting.plotting_domain(ds)
    ds_x_err = SpectralFitting.bin_widths(ds) ./ 2
    ds_data = result.objective
    ds_data_err = sqrt.(result.variance)
    ds_model = invoke_result(result, result.u)
    ds_residuals = (ds_data .- ds_model) ./ ds_data_err
    ds_residuals_err = ones(length(ds_residuals))
    # cosmetic averaging
    n = 0
    ds_x_tot = 0.0
    ds_x_binned = []
    ds_data_tot = 0.0
    ds_data_binned = []
    ds_data_err_tot = 0.0
    ds_data_err_binned = []
    ds_model_tot = 0.0
    ds_model_binned = []
    for i in 1:length(ds_x)
        n = n + 1
        ds_x_tot += ds_x[i]
        ds_data_tot += ds_data[i]
        ds_data_err_tot += ds_data_err[i]^2
        ds_model_tot += ds_model[i]
        if (i % bin_factor == 0) || (i == length(ds_x))
            push!(ds_x_binned, ds_x_tot / n)
            push!(ds_data_binned, ds_data_tot / n)
            push!(ds_data_err_binned, sqrt(ds_data_err_tot) / n)
            push!(ds_model_binned, ds_model_tot / n)
            n = 0
            ds_x_tot = 0.0
            ds_data_tot = 0.0
            ds_data_err_tot = 0.0
            ds_model_tot = 0.0
        end
    end
    ds_residuals_binned = (ds_data_binned .- ds_model_binned) ./ ds_data_err_binned
    ds_residuals_err_binned = ones(length(ds_residuals_binned))
    ds_ratio = ds_data_binned ./ ds_model_binned
    ds_ratio_err = ds_data_err_binned ./ ds_model_binned

    df = DataFrame(
        ds_e = ds_x_binned,
        ds_spec = ds_data_binned,
        ds_spec_pm = ds_data_err_binned,
        ds_model = ds_model_binned,
        ds_residuals = ds_residuals_binned,
        ds_residuals_pm = ds_residuals_err_binned,
        ds_ratio = ds_ratio,
        ds_ratio_pm = ds_ratio_err
    )
    CSV.write(file_path, df)
end

output_result(xmm, result[1], 10, "presentation/spin_results_xmm.csv")
output_result(nustar_a, result[2], 10, "presentation/spin_results_nustar_a.csv")
output_result(nustar_b, result[3], 10, "presentation/spin_results_nustar_b.csv")
