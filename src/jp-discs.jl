# fitting thin discs with non-Kerr metric
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

jp_model = JPLineProfile(a = FitParam(0.4), θ = FitParam(30.0), a13 = FitParam(0.0), eps3 = FitParam(2.0), E₀ = FitParam(6.35, frozen=true))
domain = collect(range(0.1, 1.6, 200))
output = invokemodel(domain, jp_model)
plot(domain[1:end-1], output)

model = PowerLaw(a = FitParam(2.0, lower_limit = 1.7, upper_limit = 2.2)) + AsConvolution(jp_model)(DeltaLine(E = FitParam(6.35)))

model.rin_1.frozen = true
model.rout_1.frozen = true

model

# just fit to XMM data
data_path = "data"
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





##### the following will have to be updated once the fitting is working #####

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

# output_result(xmm, result[1], 1, 2, 1, "presentation/jp_disc.csv")
