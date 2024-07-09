// Sesto 2024 presentation
// About 15 or so slides is ideal for the time available
// i.e., 15 minutes plus 5 minutes for questions

#import "@preview/polylux:0.3.1": *
#import themes.university: *

#show: university-theme.with(
  short-author: "Andrew Young",
  short-title: "Disc, corona, and spacetime models",
  short-date: "Sesto Workshop - July 2024",
)

// Slide 1
// Title slide
#title-slide(
  authors: ("Andrew Young and Fergus Baker"),
  title: "Disc, corona, and spacetime models",
  subtitle: "The flexibility to experiment with more realistic models",
  date: "July 2024",
  institution-name: "University of Bristol",
  logo: image("uob-logo.png", width: 60mm)
)

// Slide 2
#slide[
  == Introduction

  - Disc reflection and iron lines
  - Calculating line profiles for different geometries
  - _XMM-Newton_ and _NuSTAR_ observations of MCG--6-30-15
  - Measuring black holes spin
  - Thick discs with different Eddington fractions
  - Extended corona geometries
  - Testing the Kerr hypothesis by fitting "deformation parameters"
  - X-ray timing, light echoes, and reverberation
  - The future (_XRISM_, _Athena_, ...)
]

// Slide ?
// add one of the standard figures to the right - e.g., the NASA accretion disc artist's impression (we should probably make a more realistic version at some point)
#slide[
  == Questions...

  - What is the accretion disc geometry?
    - Disc thickness, inner radius
  - What is the corona geometry?
    - Size, shape, location
  - What is the spacetime geometry?
    - How fast is the black hole spinning?
    - Is General Relativity an accurate description?
]

// Slide 3
// include figure of reflection spectrum
#slide[
  == Reflection spectroscopy

  - Accretion disc is illuminated by X-rays
  - These X-rays are reprocessed ("reflected") by the disc
    - Photoelectrically absorbed (most likely below $~6-7$ keV)
      - Resulting in fluorescence ($1/3$ of the time for iron)
      - Ejection of an outer electron
    - Compton scattered (most likely above $~6-7$ keV)
  - Resulting reflection spectrum has a strong iron fluorescence line at 6.4 keV and "Compton hump" at $20-30$ keV
]

// Slide 4
// nice figure showing some ray tracing from an extended corona above the disc
#slide[
  == Calculating iron line profiles

  #side-by-side[
    - We use `Gradus.jl`
    - Disc illumination
    - Doppler shifts
    - Gravitational redshifts
    - Open source
  ][
    _Figure here_
  ]

]

// historic slide showing Fabian calculations and ASCA result - point out that thin disc, Sch / Kerr spacetime, power law emissivity adequate to describe data

// Slide 5
#slide[
  == Ratio of _XMM_ and _NuSTAR_ data to power law model

  // Describe the data we'll be using
  // perhaps make this figure a ratio plot with a wider aspect ratio

  #figure(
    image("powerlaw_fit.svg", width: 70%),
  )

  Broad iron K$alpha$ line $~6.4$ keV and Compton huump $~20-30$ keV
]

// Slide 6
// nice figure showing how disc image with changing spin for Kerr metric and corresponding line profile
// could replace one of the figures with an image showing contours of the ISCO as a function of a
#slide[
  == Spin

  #side-by-side(columns: (1fr, 1fr))[
    #figure(
      image("disc_a_0.svg", width: 60%),
    )
    #figure(
      image("disc_a_0_998.svg", width: 60%),
    )
  ][
    - Spin of the Black Hole (BH) determines the Innermost Stable Circular Orbit (ISCO)
    - Static BH, $a = 0$, $"ISCO" = 6 r_g$
    - Maximally spinning BH, $a = 0.998$, $"ISCO" = 1.235 r_g$
    - Higher spin $=>$ more extreme velocities, stronger gravitational redshift $=>$ broader iron line
  ]
]

// Slide 7
// Note that the features of the broad line are indicative of the sorts of things that can constrain models of discs, spacetimes, etc., especially when combined with NuSTAR data to constrain the continuum and Compton hump
#slide[
  == Results

  #figure(
    image("spin_results.svg", width: 60%),
  )
  Representative fit to _XMM-Newton_ data, $a = 0.998$ and $a = 0$ models.
]

// Slide 8
#slide[
  == Thick discs

  - Traditionally model a razor thin discs defined by $theta = pi/2$ as an approximation to a realistic disc with some finite thickness
  - We can model a Shakura & Sunyaev (1973) disc of finite thickness, including the relativistic effects (Page & Thorne, 1974)
  - Extends work of Taylor & Reynolds (2018)
  - In principle we can have arbitrary spacetime and corona geometries (Abdikamalov et al. 2020; see later)
  - As a proof of concept we fit one _XMM_ dataset with a thick disc of Eddington fraction $eta$ self-consistently illuminated by a lamp post corona at height $h$ around a BH of spin $a$.

  // include figure of edge-on plot showing change in thickness with Eddington fractions
  // some nice renders of disc images showing the eclipsing of the inner parts of the disc at it thickens
  // ideally would, of course, fit all data
]

// another slide showing cross sections and self-consistent lamp post illumination model
#slide[
  == Thick disc model

  #figure(
    image("disc_profile_a_0.svg", width: 60%),
  )
  Example cross section through BH and disc for $a = 0$. Lamppost is at arbitrary height. Disc thickness increases with Eddington fraction.
  // and decreases with spin
]

#slide[
  == Thick disc images
  
  Thicker discs can obscure their inner regions
  // sledging down the potential well
]

// Slide 9
// include a contour plot if possible
// also include fit parameters and error bars
#slide[
  == Results

  #figure(
    image("thick_disc.svg", width: 50%),
  )

  Thick disc fit. By eye looks the same as thin disk! However, combining _XMM_ and _NuSTAR_ data will provide stronger constraints.
]

#slide[
  == Preliminary results

  #side-by-side(columns: (1fr, 1fr))[
    #figure(
      image("thick_disc_h_eta_contours.svg", width: 100%),
      // note the degeneracy - a higher source height requires higher Eddington fraction to intercept more photons closer to the horizon
    )
  ][
    #figure(
      image("thick_disc_h_a_contours.svg", width: 100%)
      // high spin seems robustly required with little degeneracy; high source height ruled out because significant illumination of very innermost regions is required; see vertical "tail" in contours above
    )
  ]

  This simple model is consistent with low source height, $h approx 2 r_g$ few per cent of Eddington, $eta approx 4%$, high spin, $a approx 0.94$.
]


// Slide 10
#slide[
  == Corona parameters

  - Describe the model
  - Emissivity
]

// Slide 11
#slide[
  == Results

  Show the results

  // would be nice to show a plot constraining, e.g., corona height and BH spin, or similar
  // maybe corona size (w x h) and height but that might have to wait until later
]

// Slide 12
#slide[
  == Deformation parameters

  - Models typically restricted to Schwarzschild or Kerr for simplicity
  - Modified metrics have been used but required significant analytic work to determine geodesics, ISCO, redshifts, work for limited disc and corona geometries, and limited parameter flexibility
  - `Gradus.jl` allows us to work with any metric and geometry
  - Consider the "deformation" parameters that quantify departures from the Kerr metric (Johannsen 2013; Johannsen & Psaltis 2010)
  - We don't need any analytic calculations - we just feed in the metric
  - Can allow arbitrary combinations of deformation parameters not limited to varying one at a time
]

// Slide 13
#slide[
  == Results

  - We can build on the work of, e.g., Tripathi et al. (2019, 2021), Abdikamalov et al. (2020), Jiang et al. (2022)

  // would be nice to show some confidence contours of two deformation parameters :)
]

// Slide 14
// #polylux-slide[
//   == Reverberation

//   - Perhaps for the future?
// ]

// Slide 15
#slide[
  == Conclusions

  - We are now able to fit data with more realistic disc models
    - Arbitrary disc geometry
    - Arbitrary corona geometry
    - Arbitrary spacetime geometry
  - Can convolve with realistic disc reflection spectra
    - `reflionx`, `xillver`, or any other model
  - We can also do this with timing analysis
    - Not _quite_ ready in time for this talk
  - We are very happy to collaborate if you are working on anything where you think this could be helpful
]
