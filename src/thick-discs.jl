# fitting thick discs with lamp post model
# MCG--6-30-15 one XMM dataset only for illustrative purposes
# based on the spin.jl fits

using SpectralFitting
using Plots
using Gradus
using GradusSpectralModels
using JLD2
# for JLD2 compression
using CodecZlib
# using Reflionx
using LaTeXStrings
using CSV
using DataFrames

include("additional-models.jl")

function read_or_make_transfer_table(path; new=false)
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
        jldopen(path, "w"; compress=true) do f
            f["table"] = table
        end
        table
    end
end

data_path = "data"
# TRANSFER_FUNCTION_PATH = joinpath(data_path, "thin-disc-transfer-table-900.jld2")
TRANSFER_FUNCTION_PATH = joinpath(data_path, "thick-disc-transfer-table-new.jld2")

LAMP_POST_PROFILE_PATH = joinpath(data_path, "lamp-post-profile.jld2")

table = read_or_make_transfer_table(TRANSFER_FUNCTION_PATH);

for grid in table.grids
    replace!(grid.lower_f, NaN => 0.0)
    replace!(grid.upper_f, NaN => 0.0)
end

# emissivity profile table
a_range = collect(range(0.0, 0.998, 10))
h_range = collect(range(1.0, 10.0, 10))
η_range = collect(range(0.01, 0.3, 10))

prof_wrap = if isfile(LAMP_POST_PROFILE_PATH)
    f = jldopen(LAMP_POST_PROFILE_PATH, "r")
    prof = f["lamp_post"]
    close(f)
    prof
else
    prof_wrap = @time tabulate_emissivity_profile(a_range, h_range, η_range)
    jldopen(LAMP_POST_PROFILE_PATH, "w"; compress=true) do f
        f["lamp_post"] = prof_wrap
    end
end ;

# Set up model
# Note there are other options - see PR #3
lp_model = LampPostThickDisc(
    prof_wrap,
    table;
    K=FitParam(1e-4),
    θ=FitParam(30.0, lower_limit=5, upper_limit=85),
    E₀=FitParam(6.35),
    a=FitParam(0.998, lower_limit=0.0, upper_limit=0.998),
    h=FitParam(7.0, lower_limit=1.0, upper_limit=10.0),
    η=FitParam(0.15),
    # speed up model evaluation about 3 times
    quadrature_points=13,
    n_radii=600,
)

model = DeltaLine(K=FitParam(1e-5), E=FitParam(6.35, frozen=true)) + PowerLaw(K=FitParam(1e-2)) + lp_model
# model = PowerLaw(K=FitParam(1e-2)) + lp_model

domain = collect(range(1, 50, 500))
output = invokemodel(domain, model)
plot(domain[1:end-1], output, xscale=:log10, yscale=:log10)

# always use the ISCO
model.rin_1.frozen = true
model.rout_1.frozen = true
model.E₀_1.frozen = true

# try for fixed high spin
# model.a_1.frozen = true

model

# just fit to XMM data
spec = joinpath(data_path, "PN_spectrum_grp_0693781301_S003_total.fits")
back = joinpath(data_path, "PNbackground_spectrum_0693781301_S003_total.fits")
rmf = joinpath(data_path, "PN_0693781301_S003_total.rmf")
arf = joinpath(data_path, "PN_0693781301_S003_total.arf")
xmm = XmmData(XmmEPIC(), spec, background=back, response=rmf, ancillary=arf)
regroup!(xmm)
drop_bad_channels!(xmm)
mask_energies!(xmm, 3.0, 10.0)
subtract_background!(xmm)
normalize!(xmm)

# define problem and fit
prob = FittingProblem(model, xmm)

result = @time fit(prob, LevenbergMarquadt(); verbose=true, x_tol=1e-5, max_iter=100)

# plotresult(
#     xmm,
#     result[1],
#     xscale=:log10,
#     yscale=:log10,
#     xlims=(3, 30),
#     ylims=(5e-5, 2),
#     residual_ylims=(-5, 5)
# )

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
        ds_e=ds_x_binned,
        ds_spec=ds_data_binned,
        ds_spec_pm=ds_data_err_binned,
        ds_model=ds_model_binned,
        ds_residuals=ds_residuals_binned,
        ds_residuals_pm=ds_residuals_err_binned,
        ds_ratio=ds_ratio,
        ds_ratio_pm=ds_ratio_err,
        ds_model_ratio=ds_model_ratio,
        ds_no_spin_ratio=ds_no_spin_ratio
    )
    custom_header = "ds_e,ds_spec,+-,ds_model,ds_residuals,+-,ds_ratio,+-,model_ratio,no_spin"
    open(file_path, "w") do f
        println(f, custom_header)
        CSV.write(f, df, append=true, header=false)
    end
end

output_result(xmm, result[1], 1, 2, 1, "presentation/thick_disc.csv")

# create a contour plot of Eddington fraction versus source height
best_stat = result.χ2
best_h = result.u[4]
best_η = result.u[5]
h_values = collect(1.1:0.2:3.5)
η_values = collect(0.01:0.025:0.3)
Δχ² = zeros(Float64, length(η_values), length(h_values))

for (i, h) in enumerate(h_values)
    Threads.@threads for j in 1:length(η_values)
        η = η_values[j]
        _model = deepcopy(model)

        println("Fitting h = $h, η = $η")
        _model.h_1.value = h
        _model.h_1.frozen = true
        _model.η_1.value = η
        _model.η_1.frozen = true

        _prob = FittingProblem(_model, xmm)

        result = fit(_prob, LevenbergMarquadt(), x_tol=1e-4, max_iter=100)
        Δχ²[j, i] = sum(result.χ2) - best_stat
        println("Δχ² = $(Δχ²[j, i])")
    end
end

# contour levels for two free parameters
# 1, 2, and 3σ  corresponds to 68%, 95%, and 99.7% confidence intervals
# contour_levels = [2.30, 6.18, 11.83]
# corresponds to 68%, 90%, 99% confidence intervals - this is the XSPEC default
contour_levels = [0, 2.30, 4.61, 9.21]
# contour_labels = ["68%", "90%", "99%"]

contour(h_values, η_values, Δχ², levels=contour_levels, xlabel="Source height h (GM/c²)", ylabel="Eddington fraction λ", xrange=(1.1, 3.5), yrange=(0.01, 0.285), colorbar=false, color=[:green, :orange, :red], linecolor=[:black], fill=true)
scatter!([best_h], [best_η], marker=:star, markersize=16, color=:cyan, label="")
savefig("presentation/thick_disc_h_eta_contours.svg")

# create a contour plot of spin versus source height
best_stat = result.χ2
best_h = result.u[4]
best_a = result.u[2]
h_values = collect(1.1:0.2:3.5)
a_values = collect(0.9:0.01:1.0)
h_a_Δχ² = zeros(Float64, length(a_values), length(h_values))

for (i, h) in enumerate(h_values)
    Threads.@threads for j in 1:length(a_values)
        a = a_values[j]
        _model = deepcopy(model)

        println("Fitting h = $h, a = $a")
        _model.h_1.value = h
        _model.h_1.frozen = true
        _model.a_1.value = a
        _model.a_1.frozen = true

        _prob = FittingProblem(_model, xmm)

        result = fit(_prob, LevenbergMarquadt(), x_tol=1e-4, max_iter=100)
        h_a_Δχ²[j, i] = sum(result.χ2) - best_stat
        println("Δχ² = $(h_a_Δχ²[j, i])")
    end
end

# contour levels for two free parameters
# 1, 2, and 3σ  corresponds to 68%, 95%, and 99.7% confidence intervals
# contour_levels = [2.30, 6.18, 11.83]
# corresponds to 68%, 90%, 99% confidence intervals - this is the XSPEC default
contour_levels = [0, 2.30, 4.61, 9.21]
# contour_labels = ["68%", "90%", "99%"]

contour(h_values, a_values, h_a_Δχ², levels=contour_levels, xlabel="Source height h (GM/c²)", ylabel="Black hole spin a", xrange=(1.1, 3.5), yrange=(0.9, 1), colorbar=false, color=[:green, :orange, :red], linecolor=[:black], fill=true)
scatter!([best_h], [best_a], marker=:star, markersize=16, color=:cyan, label="")
savefig("presentation/thick_disc_h_a_contours.svg")
