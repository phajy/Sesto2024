# A quick look at the MCG-6-30-15 spectroscopic data

using SpectralFitting
using Plots

data_path = "data"

# XMM data

function prepare_xmm(spec, back, rmf, arf)
    ds = XmmData(
        XmmEPIC(),
        spec,
        background=back,
        response=rmf,
        ancillary=arf
    )
    regroup!(ds)
    drop_bad_channels!(ds)
    mask_energies!(ds, 3.0, 10.0)
    subtract_background!(ds)
    normalize!(ds)
    return ds
end

xmm_data = []
for i in 2:4
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
for obsid in ["60001047002"]
    for fpm in ["A", "B"]
        push!(nustar_data, prepare_nustar(data_path, obsid, fpm))
    end
end

# Fit

model = PowerLaw()
prob = FittingProblem(
    FittableMultiModel(model, model, model, model, model),
    FittableMultiDataset(xmm_data[1], xmm_data[2], xmm_data[3], nustar_data[1], nustar_data[2])
)
result = fit(prob, LevenbergMarquadt())

# Plot

plotresult(xmm_data[1], result[1], xscale=:log10, yscale=:log10, xlims=(3, 10), ylims=(5e-2, 2))
plotresult!(xmm_data[2], result[2], xscale=:log10, yscale=:log10, xlims=(3, 10), ylims=(5e-2, 2))
plotresult!(xmm_data[3], result[3], xscale=:log10, yscale=:log10, xlims=(3, 10), ylims=(5e-2, 2))
plotresult!(nustar_data[1], result[4], xscale=:log10, yscale=:log10, xlims=(3,50), ylims=(1e-4,2))
plotresult!(nustar_data[2], result[5], xscale=:log10, yscale=:log10, xlims=(3,50), ylims=(1e-4,2))
