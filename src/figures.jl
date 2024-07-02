# Figure for the presentation

using Gradus
using Plots, LaTeXStrings

# Spin

# High spin
m = KerrMetric(M=1.0, a=0.998)
x = SVector(0.0, 1000.0, deg2rad(75), 0.0)
d = ThinDisc(1.0, 20.0)

pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

α, β, img = rendergeodesics(
    m,
    x,
    d,
    # maximum integration time
    2000.0,
    αlims = (-25, 25), 
    βlims = (-12, 15),
    image_width = 800,
    image_height = 400,
    verbose = true,
    pf = pf,
)

p = heatmap(α, β, xlabel=L"\alpha \quad (r_g)", ylabel=L"\beta \quad (r_g)", img, aspect_ratio = 1)
savefig(p, "presentation/disc_a_0_998.svg")

# Low spin
m = KerrMetric(M=1.0, a=0.0)
x = SVector(0.0, 1000.0, deg2rad(75), 0.0)
d = ThinDisc(1.0, 20.0)

pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

α, β, img = rendergeodesics(
    m,
    x,
    d,
    # maximum integration time
    2000.0,
    αlims = (-25, 25), 
    βlims = (-12, 15),
    image_width = 800,
    image_height = 400,
    verbose = true,
    pf = pf,
)

p = heatmap(α, β, xlabel=L"\alpha \quad (r_g)", ylabel=L"\beta \quad (r_g)", img, aspect_ratio = 1)
savefig(p, "presentation/disc_a_0.svg")
