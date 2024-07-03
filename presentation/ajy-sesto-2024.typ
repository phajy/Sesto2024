// Get Polylux from the official package repository
#import "@preview/polylux:0.3.1": *

// Make the paper dimensions fit for a presentation and the text larger
#set page(paper: "presentation-16-9")
#set text(size: 25pt)

// About 15 or so slides is ideal for the time available
// i.e., 15 minutes plus 5 minutes for questions

// Slide 1
#polylux-slide[
  #align(horizon + center)[
    = Fitting realistic disc, corona, and \ spacetime models

    Dr Andrew Young \ University of Bristol

    With thanks to Fergus Baker, \ Jiachen Jiang & Wiktoria Tarnopolska
  ]
]

// Slide 2
#polylux-slide[
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
#polylux-slide[
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
#polylux-slide[
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
#polylux-slide[
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

// Slide 5
#polylux-slide[
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
#polylux-slide[
  == Spin

  #side-by-side[
    #figure(
      image("disc_a_0.svg", width: 60%),
      // caption: [$a = 0$ static]
    )
    #figure(
      image("disc_a_0_998.svg", width: 60%),
      // caption: [$a = 0.998$ maximally spinning]
    )
  ][
    - The spin of the Black Hole (BH) determines the Innermost Stable Circular Orbit (ISCO)
    - Static BH with $a = 0$ has ISCO at $6 r_g$
    - Maximally spinning BH with $a = 0.998$, has ISCO at $1.235 r_g$
    - More extreme velocities and stronger gravitational redshift experienced at higher spins
  ]
]

// Slide 7
#polylux-slide[
  == Results

  Show the results
]

// Slide 8
#polylux-slide[
  == Thick discs

  - Traditionally model a razor thin discs defined by $theta = pi/2$ as an approximation to a realistic disc with some finite thickness
  - We can model a Shakura & Sunyaev (1973) disc with a finite thickness, including the relativistic effects considered by Page & Thorne (1974)
  - Extends work of Taylor & Reynolds (2018)
  - In principle we can have arbitrary spacetime and corona geometries (see later)

  // include figure of edge-on plot showing change in thickness with Eddington fractions
  // some nice renders of disc images showing the eclipsing of the inner parts of the disc at it thickens
]

// Slide 9
#polylux-slide[
  == Results

  Show the results
]

// Slide 10
#polylux-slide[
  == Corona parameters

  - Describe the model
  - Emissivity
]

// Slide 11
#polylux-slide[
  == Results

  Show the results

  // would be nice to show a plot constraining, e.g., corona height and BH spin, or similar
  // maybe corona size (w x h) and height but that might have to wait until later
]

// Slide 12
#polylux-slide[
  == Deformation parameters

  - Models typically restricted to Schwarzschild or Kerr for simplicity
  - Modified metrics have been used but required significant analytic work to determine geodesics, ISCO, redshifts, work for limited disc and corona geometries, and limited parameter flexibility
  - `Gradus.jl` allows us to work with any metric and geometry
  - Consider the "deformation" parameters that quantify departures from the Kerr metric (Johannsen 2013; Johannsen & Psaltis 2010)
  - We don't need any analytic calculations - we just feed in the metric
  - Can allow arbitrary combinations of deformation parameters not limited to varying one at a time
]

// Slide 13
#polylux-slide[
  == Results

  - We can build on the work of, e.g., Tripathi et al. (2019, 2021), Abdikamalov et al. (2020), Jiang et al. (2022)

  // would be nice to show some confidence contours of two deformation parameters :)
]

// Slide 14
#polylux-slide[
  == Reverberation

  - Perhaps for the future?
]

// Slide 15
#polylux-slide[
  == Conclusions

  - Summarise the results
  - Future work
]
