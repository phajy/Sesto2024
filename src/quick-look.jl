# A quick look at the MCG-6-30-15 spectroscopic data

using SpectralFitting
using Plots

data_directory = "data"
cd(data_directory)

# load and prepare XMM obsid 0693781201
xmm_data_a = XmmData(
    XmmEPIC(),
    "PN_spectrum_grp_0693781201_S003_total.fits",
    response="PN_0693781201_S003_total.rmf",
    ancillary="PN_0693781201_S003_total.arf",
    background="PNbackground_spectrum_0693781201_S003_total.fits"
)

regroup!(xmm_data_a)
drop_bad_channels!(xmm_data_a)
mask_energies!(xmm_data_a, 3.0, 10.0)
subtract_background!(xmm_data_a)
normalize!(xmm_data_a)

# load and prepare XMM obsid 0693781301
xmm_data_b = XmmData(
    XmmEPIC(),
    "PN_spectrum_grp_0693781301_S003_total.fits",
    response="PN_0693781301_S003_total.rmf",
    ancillary="PN_0693781301_S003_total.arf",
    background="PNbackground_spectrum_0693781301_S003_total.fits"
)

regroup!(xmm_data_b)
drop_bad_channels!(xmm_data_b)
mask_energies!(xmm_data_b, 3.0, 10.0)
subtract_background!(xmm_data_b)
normalize!(xmm_data_b)

# load and prepare XMM obsid 0693781401
xmm_data_c = XmmData(
    XmmEPIC(),
    "PN_spectrum_grp_0693781401_S003_total.fits",
    response="PN_0693781401_S003_total.rmf",
    ancillary="PN_0693781401_S003_total.arf",
    background="PNbackground_spectrum_0693781401_S003_total.fits"
)

regroup!(xmm_data_c)
drop_bad_channels!(xmm_data_c)
mask_energies!(xmm_data_c, 3.0, 10.0)
subtract_background!(xmm_data_c)
normalize!(xmm_data_c)

# load and prepare NuSTAR obsid 60001047002
nustar_a_fpma = NuStarData(
    "nu60001047002A01_sr_min20.pha"
)

regroup!(nustar_a_fpma)
drop_bad_channels!(nustar_a_fpma)
mask_energies!(nustar_a_fpma, 3.0, 50.0)
subtract_background!(nustar_a_fpma)
normalize!(nustar_a_fpma)

nustar_a_fpmb = NuStarData(
    "nu60001047002B01_sr_min20.pha"
)

regroup!(nustar_a_fpmb)
drop_bad_channels!(nustar_a_fpmb)
mask_energies!(nustar_a_fpmb, 3.0, 50.0)
subtract_background!(nustar_a_fpmb)
normalize!(nustar_a_fpmb)

model = PowerLaw()
prob = FittingProblem(
    FittableMultiModel(model, model, model, model, model),
    FittableMultiDataset(xmm_data_a, xmm_data_b, xmm_data_c, nustar_a_fpma, nustar_a_fpmb)
)

result = fit(prob, LevenbergMarquadt())

plotresult(xmm_data_a, result[1], xscale=:log10, yscale=:log10, xlims=(3, 50), ylims=(1e-4, 2))
plotresult!(xmm_data_b, result[2], xscale=:log10, yscale=:log10, xlims=(3, 50), ylims=(1e-4, 2))
plotresult!(xmm_data_c, result[3], xscale=:log10, yscale=:log10, xlims=(3, 50), ylims=(1e-4, 2))
plotresult!(nustar_a_fpma, result[4], xscale=:log10, yscale=:log10, xlims=(3, 50), ylims=(1e-4, 2))
plotresult!(nustar_a_fpmb, result[5], xscale=:log10, yscale=:log10, xlims=(3, 50), ylims=(1e-4, 2))
