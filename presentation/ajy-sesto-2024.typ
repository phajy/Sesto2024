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
// #slide[
//   == Overview

//   - Disc reflection and iron lines
//   - Calculating line profiles for different geometries
//   - _XMM-Newton_ and _NuSTAR_ observations of MCG--6-30-15
//   - Measuring black holes spin
//   - Thick discs with different Eddington fractions
//   - Extended corona geometries
//   - Testing the Kerr hypothesis by fitting "deformation parameters"
//   - X-ray timing, light echoes, and reverberation
//   - The future (_XRISM_, _Athena_, ...)
// ]

// Slide ?
// add one of the standard figures to the right - e.g., the NASA accretion disc artist's impression (we should probably make a more realistic version at some point)
#slide[
  == Questions we've heard at this meeting

  - What is the accretion disc geometry?
    - Disc thickness, inner radius
  - What is the corona geometry?
    - Size, shape, location, velocity structure
  - What is the spacetime geometry?
    - How fast is the black hole spinning?
    - Is General Relativity an accurate description of spacetime?
  - Better data need better models
    - Want fast, flexible models to explore a wide range of possibilities
]

// accretion overview / history
// fabian 1988 and then asca plot (?)

#matrix-slide(columns: (2fr, 1fr))[
  == X-ray fluorescence from the inner disc

  #figure(
    image("fabian1989.png", width: 100%),
  )
  Fabian et al. (1989) $<- 35$ years ago
][
  #figure(
    image("accretion_disk.jpg", height: 90%),
  )
  #text(size: 15pt)[Credit: NASA-JPL/Caltech]
]

#slide[
  == Disc and corona

  #side-by-side(columns: (4fr, 3fr))[
    #figure(
      image("disc_and_reflection.jpg", height: 80%),
    )
  ][
    - Thermal emission from disk
    - Corona inverse-Compton scatters opt/UV $->$ X-rays
    - Disk illuminated by X-rays
    - Produces a back-scattered "reflection" spectrum
    - Spectrum modified by relativistic effects
  ]
]

// accretion overview / history
// fabian 1988 and then asca plot (?)

// Slide 3
// include figure of reflection spectrum
#slide[
  == "Reflection" spectrum in the disc's rest-frame

  #side-by-side(columns: (1fr, 1fr))[
    - Disc illuminated by X-rays
    - Absorption $->$ fluorescence
    - Scattering $->$ "Compton hump"
    - Resulting reflection spectrum has a strong iron fluorescence line at 6.4 keV and the Compton hump is at $20-30$ keV
    - Ionisation dependent (calculation is for a neutral disc)
  ][
    #figure(
      image("reflection.jpg", width: 100%),
    )
  ]
]
// include a figure of a reflection spectrum
// take from ross, fabian, young  perhaps

// Slide 4
// nice figure showing some ray tracing from an extended corona above the disc
#slide[
  == Relativistic blurring $->$ broad iron line

  #side-by-side[
    #figure(
      image("lineprofiles.ssd.png", width: 100%),
    )
    Baker & Young (in prep; Wed. talk)
  ][
    - Reflection spectrum modified by relativistic effects
    - Computed using the fast open source package `Gradus.jl`
    - Arbitrary disc, corona, spacetime geometries
    - Figure shows line profiles from thin discs with different black hole spin and disc inclinations
  ]
]

#slide[
  == Blurred reflection spectrum

  #figure(
    image("blurred_reflection.svg", height: 85%),
  )
]

// historic slide showing Fabian calculations and ASCA result - point out that thin disc, Sch / Kerr spacetime, power law emissivity adequate to describe data

// Slide 5
#slide[
  == MCG--6-30-15: _XMM_ & _NuSTAR_ ratio to power law

  #figure(
    image("powerlaw_fit.svg", width: 60%),
  )

  Want to fit data. Here's are some nice "clean" datasets for MCG-6-30-15. Broad iron K$alpha$ line $~6.4$ keV and Compton hump $~20-30$ keV.
]

// Slide 6
// nice figure showing how disc image with changing spin for Kerr metric and corresponding line profile
// could replace one of the figures with an image showing contours of the ISCO as a function of a
#slide[
  == Model spinning BH with a razor thin disc

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
  == Results -- ratio to continuum plus distant Fe

  #figure(
    image("spin_results.svg", width: 55%),
  )
  Representative fit to _XMM-Newton_ data, $a = 0.998$ and $a = 0$ models. High spin required to fit the extremely broad red wing of the iron line.
]

// Slide 8
#slide[
  == Better data need better models -- thick discs

  - Traditionally model a razor thin discs defined by $theta = pi/2$ as an approximation to a realistic disc with some finite thickness
  - We can model a (Shakura & Sunyaev 1973; Novikov & Thorne 1973; Page & Thorne, 1974) with different accretion rates, $dot(m)$
  - Extends the work of, e.g., Taylor & Reynolds (2018)
  // - In principle we can have arbitrary spacetime and corona geometries (e.g., Abdikamalov et al. 2020; see later)
  - As a proof of concept we again fit one _XMM_ dataset with a thick disc of Eddington fraction $accent(M, dot) / accent(M_"Edd", dot)$ self-consistently illuminated by a lamp post corona at height $h$ around a black hole of spin $a$.
  // - We should really fit the combined _XMM_ and _NuSTAR_ data to get stronger constraints (not ready yet; in preparation)
]

#slide[
  == Thick disc model

  #figure(
    image("disc_profile_a_0.svg", width: 55%),
  )
  Cross section through BH and disc for $a = 0$. Lamp post is at arbitrary height, $h$. Disc thickness increases with Eddington fraction. Illumination is calculated self-consistently, not assumed to be a power-law.
 ]

#slide[
  == Sledging down the potential well
  
  #figure(
    image("transfer-function.parameterization.svg", width: 100%),
  )

  A thick disc (b) is much better for sledging down the potential well than a razor thin disc (a). At higher inclinations thick discs eclipse their inner regions, including the ISCO (thick black line).
]

#slide[
  == Self-eclipsing of the disc

  #figure(
    image("datum-plane.svg", width: 60%),
  )

  More important at high inclination angles or for very thick discs (may be relevant, e.g., for tidal disruption events, super-Eddington sources).
]

// Slide 9
// include a contour plot if possible
// also include fit parameters and error bars
#slide[
  == Example thick disc fit

  #figure(
    image("thick_disc.svg", width: 50%),
  )

  By eye looks the same as thin disk! However, even with just the _XMM_ data we can still constrain the accretion flow geometry.
]

#slide[
  == Preliminary results (more data $->$ better constraints)

  #side-by-side(columns: (1fr, 1fr))[
    #figure(
      image("thick_disc_h_eta_contours.svg", width: 100%),
      // note the degeneracy - a higher source height requires higher Eddington fraction to intercept more photons closer to the horizon
    )
  ][
    #figure(
      image("thick_disc_h_a_contours.svg", width: 100%)
      // high spin seems robustly required with little degeneracy; high source height ruled out because significant illumination of very innermost regions is required; see vertical "tail" in contours above
      // high spin also has thinner disc at given eta
    )
  ]

  Low source height, $h approx 2 r_g$, Eddington fraction $approx 4%$, high spin, $a approx 0.94$. Consistent with results of Jiang et al. (2021).
]

// Slide 10
// #slide[
//   == Corona parameters

//   - Describe the model
//   - Emissivity
// ]

// Slide 11
// #slide[
//   == Results

//   Show the results

  // would be nice to show a plot constraining, e.g., corona height and BH spin, or similar
  // maybe corona size (w x h) and height but that might have to wait until later
// ]

// Slide 12
#slide[
  == Deformation parameters -- non-standard spacetimes

  - Models typically restricted to Schwarzschild or Kerr for simplicity
  - Modified metrics have been used but required significant analytic work to determine geodesics, ISCO, redshifts, work for limited disc and corona geometries, and limited parameter flexibility
  - `Gradus.jl` allows us to work with any metric and geometry
  - Consider the "deformation" parameters that quantify departures from the Kerr metric (Johannsen 2013; Johannsen & Psaltis 2010)
  - We don't need any analytic calculations -- we just feed in the metric
  - Can allow arbitrary combinations of deformation parameters not limited to varying one at a time
]

// Slide 13
#slide[
  == Results

  - We can build on the work of, e.g., Tripathi et al. (2019, 2021), Abdikamalov et al. (2020), Jiang et al. (2022)
  - Can have thick disc, self-consistently illuminated, with any spacetime
  - Will extend to arbitrary corona geometry

  // would be nice to show some confidence contours of two deformation parameters :)
  // could note that a major constraint here is the isco which is deformation parameter dependent too
  // note the tarnopolska paper in prep about self-consistent thick discs with deformation parameters
]

// describe how line profiles change with deformation parameters

// Slide 14
// #polylux-slide[
//   == Reverberation

//   - Perhaps for the future?
// ]

// include a real-life demo at the end of the talk? perhaps compute a line profile for complex case convolved with reflionx and show the result, maybe even fit to the data ... ?
// note that timing properties are important for breaking these degeneracies - light travel time to narrower components for example - thin disc longer lag than thick discs

// include cosimo summary plot of deformation constraints? or compare with eht which is in weaker field - not probing so close to the horizon which is critically important for being so ultra-deep in the metric (photon ring not deep enough) - unique to x-rays within 2 r_g (!) - state photon rin radius as function of a - plot of relevant radii; show gradus photon ring plot; test in *strong* field regime

// Slide 15
#slide[
  == Conclusions

  - We are now able to fit data with more realistic disc models
    - Arbitrary disc geometry
    - Arbitrary corona geometry
    - Arbitrary spacetime geometry
  - Can convolve with realistic disc reflection spectra
    - `reflionx`, `xillver`, or any other model, including high density
  - TODO: timing studies (lags versus frequency, energy)
  - TODO: radiative transfer -- absorption, scattering, and polarisation
  - We are *very happy to collaborate* if this could be useful, or you want a comparison with another code; all our software is open source
  // although work in progress so probably easiest to talk to us, but you can raise an issue on GitHub and we'll respond
]
