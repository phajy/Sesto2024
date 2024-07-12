# fitting thin discs with non-Kerr metric
# MCG--6-30-15 one XMM dataset only for illustrative purposes

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

function read_transfer_path(path)
    if isfile(path)
        @info "Reading transfer function table from file..."
        f = jldopen(path, "r")
        table = f["table"]
        close(f)

        N = 10
        a_range = collect(range(0.0, 0.998, N))
        θ_range = collect(range(20.0, 60.0, N))
        eps_range = collect(range(-0.5, 3, N))
        Gradus.CunninghamTransferTable((a_range, θ_range, eps_range), table)
    else
        @error "Not a file"
    end
end

data_path = "data"
table = read_transfer_path(joinpath(data_path, "johannsen_transfer_table_extended.jld2"));

for grid in table.grids
    replace!(grid.lower_f, NaN => 0.0)
    replace!(grid.upper_f, NaN => 0.0)
end

profile = begin
    f = jldopen(joinpath(data_path, "johannsen_lamp_post_extended.jld2"), "r")
    prof = f["lamp_post"]
    close(f)
    prof
end;

# jp_model = JohannsenPsaltisLampPost(profile, table)
jp_model = JohannsenPsaltisFixed(table)

jp_model.K.value = 1e-2
jp_model.a.value = 0.2
# jp_model.eps.value = 0.0
jp_model.eps.lower_limit = -0.5
jp_model.eps.upper_limit = 3.0
jp_model.E₀.value = 6.35
# jp_model.Γ.value = 3.0
# jp_model.h.value = 10.0
# jp_model.h.lower_limit = 1.0
jp_model.rout.value = 300.0
# jp_model.eps.frozen = false
jp_model.rin.frozen = true
jp_model.rout.frozen = true
jp_model.E₀.frozen = true
# jp_model.Γ.frozen = false

domain = collect(range(1, 10, 200))
output = invokemodel(domain, jp_model)

plot(domain[1:end-1], output)

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
model = DeltaLine(K=FitParam(1e-4), E=FitParam(6.35, frozen=true)) + PowerLaw(K=FitParam(0.1), a = FitParam(2.0, lower_limit = 1.7, upper_limit = 2.2)) + jp_model

output = invokemodel(domain, model)
plot(domain[1:end-1], output)

prob = FittingProblem(model, xmm)

result = @time fit!(prob, LevenbergMarquadt(); verbose=true, x_tol=1e-7, max_iter=100)

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

output_result(xmm, result[1], 1, 4, 1, "presentation/jp_disc.csv")

# create a contour plot 
best_stat = result.χ2
best_a = result.u[2]
best_ϵ = result.u[4]
a_values = collect(0.5:0.1:1.0)
ϵ_values = collect(-0.5:0.2:3.0)
Δχ² = zeros(Float64, length(ϵ_values), length(a_values))

# now we've fit it, freeze the emissivity index
# model.Γ_1.frozen = true

for (i, a) in enumerate(a_values)
    Threads.@threads for j in 1:length(ϵ_values)
        ϵ = ϵ_values[j]
        _model = deepcopy(model)

        println("Fitting a = $a, ϵ = $ϵ")
        _model.a_1.value = a
        _model.a_1.frozen = true
        _model.eps_1.value = ϵ
        _model.eps_1.frozen = true

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

contour(a_values, ϵ_values, Δχ², levels=contour_levels, xlabel="a", ylabel="ϵ", xrange=(0.5, 1.0), yrange=(-0.5, 2.9), colorbar=false, color=[:green, :orange, :red], linecolor=[:black], fill=true)
# scatter!([best_a], [best_ϵ], marker=:star, markersize=16, color=:cyan, label="")
scatter!([best_a], [best_ϵ], marker=:star, markersize=16, color=:cyan, label="")

# calculate exclusion regions
function calc_exclusion(as, ϵs)
    regions = zeros(Float64, (length(as), length(ϵs)))
    Threads.@threads for i in eachindex(as)
        a = as[i]
        for (j, ϵ) in enumerate(ϵs)
            m = JohannsenPsaltisMetric(M = 1.0, a = a, ϵ3 = ϵ)
            regions[i, j] = if is_naked_singularity(m)
                NaN
            else
                Gradus.isco(m)
            end
        end
    end
    regions
end

as = range(0.0, 1.0, 100)
ϵs = range(-0.5, 3.0, 100)

img = calc_exclusion(as, ϵs)

lim_ϵ = zeros(length(as))
for (i, a) in enumerate(as)
    lim_ϵ[i] = 5.0
    for (j, ϵ) in enumerate(ϵs)
        if isnan(img[i, j])
            lim_ϵ[i] = ϵ
            break
        end
    end
end

plot!(as, lim_ϵ, color=:black, linewidth=4, label="Limit", legend=true)
plot!([0.5, 1.0], [0, 0], color=:blue, linewidth=2, label="Kerr")

savefig("presentation/jp_disc_a_eps_contours.svg")
