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

    Thanks to Fergus Baker & Jiachen Jiang
    // and others
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
#polylux-slide[
  == Calculating iron line profiles

  Something about Gradus.jl
]

// Slide 5
#polylux-slide[
  == XMM and NuSTAR data

  Describe the data we'll be using
]

// Slide 6
#polylux-slide[
  == Spin

  Describe the model
]

// Slide 7
#polylux-slide[
  == Results

  Show the results
]

// Slide 8
#polylux-slide[
  == Thick discs

  Describe the model
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
]

// Slide 12
#polylux-slide[
  == Deformation parameters

  - Describe the model
  - Multiple parameters
]

// Slide 13
#polylux-slide[
  == Results

  Show the results
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
