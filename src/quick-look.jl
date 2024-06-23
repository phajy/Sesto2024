# A quick look at the MCG-6-30-15 spectroscopic data

using SpectralFitting
using Plots
using Gradus
using GradusSpectralModels
using JLD2
# for JLD2 compression
using CodecZlib

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
TRANSFER_FUNCTION_PATH = joinpath(data_path, "transfer_function_table.jld2")

table = read_or_make_transfer_table(TRANSFER_FUNCTION_PATH);

# XMM data

function prepare_xmm(spec, back, rmf, arf)
    ds = XmmData(XmmEPIC(), spec, background = back, response = rmf, ancillary = arf)
    regroup!(ds)
    drop_bad_channels!(ds)
    mask_energies!(ds, 3.0, 10.0)
    subtract_background!(ds)
    normalize!(ds)
    return ds
end

xmm_data = []
for i = 2:4
    spec = joinpath(data_path, "PN_spectrum_grp_0693781$(i)01_S003_total.fits")
    back = joinpath(data_path, "PNbackground_spectrum_0693781$(i)01_S003_total.fits")
    rmf = joinpath(data_path, "PN_0693781$(i)01_S003_total.rmf")
    arf = joinpath(data_path, "PN_0693781$(i)01_S003_total.arf")
    push!(xmm_data, prepare_xmm(spec, back, rmf, arf))
end

# NuSTAR data

function prepare_nustar(data_path, obsid, fpm)
    ds = NuStarData(joinpath(data_path, "nu$(obsid)$(fpm)01_sr_min20.pha"))
    regroup!(ds)
    drop_bad_channels!(ds)
    mask_energies!(ds, 3.0, 50.0)
    subtract_background!(ds)
    normalize!(ds)
    return ds
end

nustar_data = []
for obsid in ["60001047002", "60001047003", "60001047005"]
    for fpm in ["A", "B"]
        push!(nustar_data, prepare_nustar(data_path, obsid, fpm))
    end
end


# Fit

model =
    PowerLaw() + LineProfile(
        x -> x^-3,
        table;
        E₀ = FitParam(6.4),
        a = FitParam(0.998, lower_limit = 0.0, upper_limit = 0.998),
    )
domain = collect(range(1, 10, 200))
output = invokemodel(domain, model)

plot(domain[1:end-1], output)

# always use the ISCO
model.rin_1.frozen = true
model.rout_1.frozen = true
model.E₀_1.frozen = true

model

datasets = FittableMultiDataset(
    xmm_data[1],
    xmm_data[2],
    xmm_data[3],
    nustar_data[1],
    nustar_data[2],
    nustar_data[3],
    nustar_data[4],
    nustar_data[5],
    nustar_data[6],
)
datasets = FittableMultiDataset(xmm_data[1], xmm_data[2], xmm_data[3])

prob = FittingProblem(model => xmm_data[3])

# bind!(prob, :a_2)
# bind!(prob, :a_1)
# bind!(prob, :θ_1)
# bind!(prob, :rout_1)
# bind!(prob, :E₀_1)
# TODO: figure out how to do the bindings properly
# bind!(prob, 4 => :K, 5 => :K)

result = fit!(prob, LevenbergMarquadt(); verbose = true, autodiff = :finite)
model

# Plot
plotresult(
    xmm_data[3],
    result[1],
    xscale = :log10,
    yscale = :log10,
    xlims = (3, 10),
    ylims = (5e-2, 2),
)

plotresult(
    xmm_data[1],
    result[1],
    xscale = :log10,
    yscale = :log10,
    xlims = (3, 10),
    ylims = (5e-2, 2),
)
plotresult!(
    xmm_data[2],
    result[2],
    xscale = :log10,
    yscale = :log10,
    xlims = (3, 10),
    ylims = (5e-2, 2),
)
plotresult!(
    xmm_data[3],
    result[3],
    xscale = :log10,
    yscale = :log10,
    xlims = (3, 10),
    ylims = (5e-2, 2),
)
plotresult!(
    nustar_data[1],
    result[4],
    xscale = :log10,
    yscale = :log10,
    xlims = (3, 50),
    ylims = (1e-4, 2),
)
plotresult!(
    nustar_data[2],
    result[5],
    xscale = :log10,
    yscale = :log10,
    xlims = (3, 50),
    ylims = (1e-4, 2),
)
plotresult!(
    nustar_data[3],
    result[6],
    xscale = :log10,
    yscale = :log10,
    xlims = (3, 50),
    ylims = (1e-4, 2),
)
plotresult!(
    nustar_data[4],
    result[7],
    xscale = :log10,
    yscale = :log10,
    xlims = (3, 50),
    ylims = (1e-4, 2),
)
plotresult!(
    nustar_data[5],
    result[8],
    xscale = :log10,
    yscale = :log10,
    xlims = (3, 50),
    ylims = (1e-4, 2),
)
plotresult!(
    nustar_data[6],
    result[9],
    xscale = :log10,
    yscale = :log10,
    xlims = (3, 50),
    ylims = (1e-4, 2),
)
