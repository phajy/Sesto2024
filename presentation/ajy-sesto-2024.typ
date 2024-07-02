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

// Slide 3
#polylux-slide[
  == Iron lines

  How are iron lines produced
]

// Slide 4
// nice figure showing some ray tracing from an extended corona above the disc
#polylux-slide[
  == Calculating iron line profiles

  #side-by-side[
    - Something about Gradus.jl
    - Open source
  ][
    _Figure here_
  ]

]

// Slide 5
#polylux-slide[
  == _XMM_ and _NuSTAR_ data

  Describe the data we'll be using
]

// Slide 6
// nice figure showing how disc image with changing spin for Kerr metric and corresponding line profile
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
