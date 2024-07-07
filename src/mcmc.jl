# MCMC test
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
using CSV
using DataFrames

using StatsPlots
using Turing

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
end

# Set up model
# Note there are other options - see PR #3
lp_model = LampPostThickDisc(
    prof_wrap,
    table;
    K=FitParam(1e-4),
    θ=FitParam(30.0, lower_limit=5, upper_limit=85),
    E₀=FitParam(6.35),
    a=FitParam(0.998, lower_limit=0.0, upper_limit=0.998),
    h=FitParam(4.0, lower_limit=1.0, upper_limit=10.0),
    # speed up model evaluation about 3 times
    quadrature_points=13,
    n_radii=600,
)

model = DeltaLine(K=FitParam(1e-5), E=FitParam(6.35, frozen=true)) + PowerLaw(K=FitParam(1e-2)) + AutoCache(lp_model)

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

result = @time SpectralFitting.fit(prob, LevenbergMarquadt(); verbose=true, x_tol=1e-5, max_iter=100)

plotresult(
    xmm,
    result[1],
    xscale=:log10,
    yscale=:log10,
    xlims=(3, 10),
    ylims=(5e-5, 2),
    residual_ylims=(-5, 5) # <-- doesn't work as expected
)

# set up MCMC fit
@model function mcmc_model(domain, objective, variance, f)
    K_1 ~ truncated(Normal(1e-5, 5e-6); lower = 0)
    a_1 ~ truncated(Normal(0.998, 0.1); lower = 0, upper = 1)
    θ_1 ~ Normal(35, 10)
    h_1 ~ truncated(Normal(4.0, 0.4); lower = 1, upper = 10)
    η_1 ~ truncated(Normal(0.1, 0.05); lower = 0, upper = 1)
    K_2 ~ truncated(Normal(0.01, 0.005); lower = 0)
    a_2 ~ Normal(2.0, 0.2)
    K_3 ~ truncated(Normal(1e-5, 5e-6); lower = 0)

    pred = f(domain, [K_1, a_1, θ_1, h_1, η_1, K_2, a_2, K_3])
    return objective ~ MvNormal(pred, sqrt.(variance))
end

config = FittingConfig(FittingProblem(model => xmm))
mm = mcmc_model(
    make_model_domain(ContiguouslyBinned(), xmm),
    make_objective(ContiguouslyBinned(), xmm),
    make_objective_variance(ContiguouslyBinned(), xmm),
    # _f_objective returns a function used to evaluate and fold the model through the data
    SpectralFitting._f_objective(config),
)

chain = sample(mm, NUTS(), 5_000)

# corner plot should be compared with the contour plot generated separately
