# A quick look at the MCG-6-30-15 spectroscopic data

using SpectralFitting
using Plots

data_directory = "data"
cd(data_directory)
xmm_data_a = XmmData(XmmEPIC(), "PN_spectrum_grp_0693781201_S003_total.fits", response="PN_0693781201_S003_total.rmf", ancillary="PN_0693781201_S003_total.arf", background="PNbackground_spectrum_0693781201_S003_total.fits")

regroup!(xmm_data_a)
drop_bad_channels!(xmm_data_a)
mask_energies!(xmm_data_a, 3.0, 10.0)
subtract_background!(xmm_data_a)
normalize!(xmm_data_a)

model = PowerLaw()
prob = FittingProblem(model => xmm_data_a)
result = fit(prob, LevenbergMarquadt())

plotresult(xmm_data_a, result, xscale=:log10, yscale=:log10, xlims=(3, 10), ylims=(0.05, 2))
