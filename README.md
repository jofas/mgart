![mgart](./icon.svg)

# mgart

## TODO

* [x] `julia set`: scale

* [x] `julia set`: zoom around a point

* [x] 1d color maps

* [x] `julia set`: coloring

* [x] `julia set`: smooth coloring

* [x] unified cli instead of binaries

* [x] mandelbrot

* [x] from config file cli option

* [x] double precision

* [x] `julia set`: parallelize

* [x] `julia set` progress report

* [x] `julia set` optimize

* [x] 1d color map: sine

* [x] `julia set`: LCH coloring: LCH to RGB (LCH -> LAB -> RGB)

* [x] `julia set`: LCH coloring: LCH interpolation

* [x] `julia set`: option for color space (rgba or lch)

* [x] fix color user api

* [x] no flatten buffer but a pixel is `chunks_mut(3)`

* [x] all gradients computed in LCH space

* [x] RGBA -> RGB

* [x] `julia set`: optimize inner loop

* [x] LCH coloring: RGB to LCH (RGB -> LAB -> LCH)

* [x] reimplement with `num_complex`

* [x] instead of color map with gradient, compute color directly (LCH)

* [x] buddhabrot

* [x] NLM with global h

* [x] PostProcessing

* [ ] CLAHE

* [ ] time estimates

* [ ] API for describing drawing algorithms

* [ ] drawing algorithms: 
  
      outer: 
      - [x] EscapeTime
      - [ ] OuterDistance
    
      inner:
      - [ ] FiniteAttractors (magnitude `dz0`)
      - [x] InnerDistance

      both:
      - [ ] Constant
      - [ ] EscapeAngle
      - [ ] CurvatureEstimation
      - [ ] Atoms

      else:
      - [ ] PickoverStalks

* [ ] smarter sampling (less iterations)

* [ ] sampling visualization

* [ ] rotation parameter

* [ ] color algorithms: ColorMap1D, LCH, RGB, Nebula

* [ ] color smoothing optional

* [ ] gradients for distance instead of [0, 1]

* [ ] wasm based editor

* [ ] find cool points and colormaps where I can see the effect of
  the following visualization methods

* [ ] curvature estimation

* [ ] statistical methods (mean, avg, etc.)

* [ ] `julia set`: field lines

* [ ] progressive iteration visualization as GIF

* [ ] newton fractals

* [ ] multibrot

* [ ] fractal flame

* [ ] more gradients

* [ ] arbitrary precision types for higher zoom levels
