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

# longest XMM and NuSTAR data sets for the sake of example
# will use _all_ data sets for publication
# attempted to use merged data for illustrative purposes but that was not a good idea

# XMM data
spec = joinpath(data_path, "PN_spectrum_grp_0693781301_S003_total.fits")
back = joinpath(data_path, "PNbackground_spectrum_0693781301_S003_total.fits")
rmf = joinpath(data_path, "PN_0693781301_S003_total.rmf")
arf = joinpath(data_path, "PN_0693781301_S003_total.arf")
xmm = XmmData(XmmEPIC(), spec, background = back, response = rmf, ancillary = arf)
regroup!(xmm)
drop_bad_channels!(xmm)
mask_energies!(xmm, 3.0, 10.0)
subtract_background!(xmm)
normalize!(xmm)

# NuSTAR data
nustar_a = NuStarData(joinpath(data_path, "nu60001047003A01_sr_min20.pha"))
regroup!(nustar_a)
drop_bad_channels!(nustar_a)
mask_energies!(nustar_a, 3.0, 50.0)
subtract_background!(nustar_a)
normalize!(nustar_a)

nustar_b = NuStarData(joinpath(data_path, "nu60001047003B01_sr_min20.pha"))
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
    θ = FitParam(35.0, lower_limit = 10.0, upper_limit = 80.0),
    rout = FitParam(400.0, frozen=true),
    # speed up model evaluation about 3 times
    quadrature_points = 13,
    n_radii = 600,
)

# power law model - used in the presentation to show the iron line and Compton hump
# model = PowerLaw()

# use the AutoCache wrapper to avoid re-evaluating an expensive model unnecessary

# model = PowerLaw(a = FitParam(2.0, lower_limit = 1.7, upper_limit = 2.2)) + AsConvolution(lp_model)(ReflionxTable(K = FitParam(1e-7), refl_table))

# model = XS_CutOffPowerLaw(Γ = FitParam(2.0, lower_limit = 1.7, upper_limit = 2.2), Ecut = FitParam(125.0, lower_limit=100.0, upper_limit=200.0), z = FitParam(0.007749, frozen=true)) + AsConvolution(lp_model)(ReflionxTable(K = FitParam(1e-7), logξ = FitParam(3.0), refl_table))

# model = XS_CutOffPowerLaw(Ecut = FitParam(125.0, lower_limit=100.0, upper_limit=200.0), z = FitParam(0.007749, frozen=true))

# comparison with Laor model
# model = PowerLaw(K = FitParam(1.0E-3)) + XS_Laor(K = FitParam(1.0E-5), lineE = FitParam(6.35, frozen=true), θ = FitParam(30.0, lower_limit = 10.0, upper_limit = 80.0))

# simple fit with iron line and r^-3 emissivity profile
# we can fit this to just the xmm data set
model = DeltaLine(K = FitParam(1e-5), E = FitParam(6.35, frozen = true)) + PowerLaw(K = FitParam(1.0E-3, upper_limit=0.1), a = FitParam(2.0, lower_limit=1.7, upper_limit=2.2)) + AsConvolution(lp_model, domain = collect(range(0, 2, 500)))(DeltaLine(K = FitParam(1.0E-5, upper_limit=0.1), E = FitParam(6.35, frozen=true, lower_limit=6.0, upper_limit=7.0)))

domain = collect(range(1, 50, 500))

output = invokemodel(domain, model)
plot(domain[1:end-1], output, xscale=:log10, yscale=:log10)

# always use the ISCO
model.rin_1.frozen = true
# model.rout_1.frozen = true
model.E₀_1.frozen = true

# model.a_1.frozen = true
# model.θ_1.frozen = true

# model.logξ_1.frozen = true
# model.E_cut_1.frozen = true
# model.Γ_1.frozen = true

# make sure the datasets from the same observatory are grouped together
# else the AutoCache will trigger re-eval as the domain has changed, even through the model parameters will all be the same
# datasets = FittableMultiDataset(
#     xmm, nustar_a, nustar_b
# )
# models = FittableMultiModel((model for _ in datasets.d)...)

# prob = FittingProblem(models, datasets)

# just fit to XMM
prob = FittingProblem(model, xmm)

# cut off power law
# bind!(prob, :Γ)
# bind!(prob, :Ecut)

# power law only
# bind!(prob, :a)

# reflection model
# bind!(prob, :a_2)
# bind!(prob, :a_1)
# bind!(prob, :θ_1)
# bind!(prob, :Γ_1)
# bind!(prob, :logξ_1)
# bind!(prob, :E_cut_1)

# cutoff powerlaw reflection
# bind!(prob, :Γ_1)
# bind!(prob, :θ_1)
# bind!(prob, :Γ_2)
# bind!(prob, :Ecut_1)

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
function output_result(ds, result, zero_index, spin_index, bin_factor, file_path)
    ds_x = SpectralFitting.plotting_domain(ds)
    ds_x_err = SpectralFitting.bin_widths(ds) ./ 2
    ds_data = result.objective
    ds_data_err = sqrt.(result.variance)
    ds_model = copy(invoke_result(result, result.u))

    # remove line
    tmp = result.u[zero_index]
    result.u[zero_index] = 0.0
    ds_model_no_line = copy(invoke_result(result, result.u))
    result.u[zero_index] = tmp

    # zero spin line
    tmp = result.u[spin_index]
    result.u[spin_index] = 0.0
    ds_model_no_spin = copy(invoke_result(result, result.u))
    result.u[spin_index] = tmp

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
    ds_model_no_line_tot = 0.0
    ds_model_no_line_binned = []
    ds_model_no_spin_tot = 0.0
    ds_model_no_spin_binned = []
    for i in 1:length(ds_x)
        n = n + 1
        ds_x_tot += ds_x[i]
        ds_data_tot += ds_data[i]
        ds_data_err_tot += ds_data_err[i]^2
        ds_model_tot += ds_model[i]
        ds_model_no_line_tot += ds_model_no_line[i]
        ds_model_no_spin_tot += ds_model_no_spin[i]
        if (i % bin_factor == 0) || (i == length(ds_x))
            push!(ds_x_binned, ds_x_tot / n)
            push!(ds_data_binned, ds_data_tot / n)
            push!(ds_data_err_binned, sqrt(ds_data_err_tot) / n)
            push!(ds_model_binned, ds_model_tot / n)
            push!(ds_model_no_line_binned, ds_model_no_line_tot / n)
            push!(ds_model_no_spin_binned, ds_model_no_spin_tot / n)
            n = 0
            ds_x_tot = 0.0
            ds_data_tot = 0.0
            ds_data_err_tot = 0.0
            ds_model_tot = 0.0
            ds_model_no_line_tot = 0.0
            ds_model_no_spin_tot = 0.0
        end
    end
    ds_residuals_binned = (ds_data_binned .- ds_model_binned) ./ ds_data_err_binned
    ds_residuals_err_binned = ones(length(ds_residuals_binned))
    ds_ratio = ds_data_binned ./ ds_model_no_line_binned
    ds_ratio_err = ds_data_err_binned ./ ds_model_no_line_binned
    ds_model_ratio = ds_model_binned ./ ds_model_no_line_binned
    ds_no_spin_ratio = ds_model_no_spin_binned ./ ds_model_no_line_binned
    ds_no_spin_ratio = (ds_no_spin_ratio .- 1.0) .* (maximum(ds_model_ratio) - 1.0) / (maximum(ds_no_spin_ratio) - 1.0) .+ 1.0

    df = DataFrame(
        ds_e = ds_x_binned,
        ds_spec = ds_data_binned,
        ds_spec_pm = ds_data_err_binned,
        ds_model = ds_model_binned,
        ds_residuals = ds_residuals_binned,
        ds_residuals_pm = ds_residuals_err_binned,
        ds_ratio = ds_ratio,
        ds_ratio_pm = ds_ratio_err,
        ds_model_ratio = ds_model_ratio,
        ds_no_spin_ratio = ds_no_spin_ratio
    )
    custom_header = "ds_e,ds_spec,+-,ds_model,ds_residuals,+-,ds_ratio,+-,model_ratio,no_spin"
    open(file_path, "w") do f
        println(f, custom_header)
        CSV.write(f, df, append=true, header=false)
    end
end

output_result(xmm, result[1], 1, 2, 1, "presentation/spin_results_xmm.csv")
# output_result(nustar_a, result[2], 10, "presentation/spin_results_nustar_a.csv")
# output_result(nustar_b, result[3], 10, "presentation/spin_results_nustar_b.csv")
