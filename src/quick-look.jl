# A quick look at the MCG-6-30-15 spectroscopic data

using SpectralFitting
using Plots

data_directory = "data"
cd(data_directory)

# load and prepare obsid 0693781201
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

# load and prepare obsid 0693781301
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

# load and prepare obsid 0693781401
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

model = PowerLaw()
prob = FittingProblem(
    FittableMultiModel(model, model, model),
    FittableMultiDataset(xmm_data_a, xmm_data_b, xmm_data_c)
)

result = fit(prob, LevenbergMarquadt())

plotresult(xmm_data_a, result[1], xscale=:log10, yscale=:log10, xlims=(3, 10), ylims=(0.05, 2))
plotresult!(xmm_data_b, result[2], xscale=:log10, yscale=:log10, xlims=(3, 10), ylims=(0.05, 2))
plotresult!(xmm_data_c, result[3], xscale=:log10, yscale=:log10, xlims=(3, 10), ylims=(0.05, 2))
