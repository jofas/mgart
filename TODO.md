# TODO

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

* [x] CLAHE

* [x] buddhabrot exponent param

* [x] sampling api

* [x] uniform random sampler

* [x] weighted KDE vs KDE

* [x] fix non-local means 

* [x] CLAHE fixed north and west borders

* [x] `Uniform` sampler: `r` parameter

* [x] `post_processing` module (with submodules)

* [x] `coloring` module (colors as submodule)

* [x] `gradient` module

* [x] rotation parameter (rotate viewport, not space)

* [x] api integration test (using example files)

* [x] cli input file as positional argument

* [x] libsonnet files (commands, types, constants)

* [x] SamplerArgs into sampler module

* [x] weighted KDE: stratified sampling

* [x] cli: add jsonnet file support to input file 

* [x] command => algorithm

* [x] build action

* [x] badges

* [x] top-level basic documentation 

* [x] `Algorithms` object-oriented 

* [x] use `unreachable` macro

* [x] CI benchmarking setup

* [x] run actions with minimal permissions

* [x] codecov comment on PRs

* [x] cache rust builds 

* [x] nicer progress printing abstraction

* [x] split image generation from drawing algorithms

* [x] Frame object (with indexing and iteration)

* [x] make SAT nice (minimize casting)

* [x] Viewport translate grid_pos to point

* [x] Viewport rotated and unrotated version of each method

* [x] benchmark buddhabrot creation 

* [x] optimize buddhabrot `iter_mandel_check_vp` (rename to 
  `trace_point`)

* [x] sampling benchmark

* [x] make benchmarking more resilient (aim for 1s time with 10
  samples)

* [x] optimize buddhabrot sampling

* [x] creator pattern for each algorithm

* [ ] get code to pass pedantic clippy with as little allows as 
  possible

* [ ] replace `divrem` and `num_complex` with `num`

* [ ] increase coverage

* [ ] collect stats (hits, avg. iterations, etc.)

* [ ] kde sampling: make sure sampling points are evenly distributed
  in space 

* [ ] api integration test for jsonnet files 
  (using example files - `jsonnet-rs v0.17.0`)

* [ ] buddhabrot: flag for disabling `iter_mande_check_vp` (makes more
  restrictive samplers faster)

* [ ] examples (jsonnet, zoom, samplers, colored)

* [ ] publish action (publish on crates.io + create zip 
  with libsonnet files)

* [ ] publish `0.0.2`

* [ ] anti-buddhabrot (escape cycles using `period` and 
  `finite_attractor`, see interior distance rendering of julia set)

* [ ] color algorithms: ColorMap1D, LCH, RGB, Nebula

* [ ] escape function parameter (try Rhai)

* [ ] clean-up julia-set implementation

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

* [ ] multibrot

* [ ] sampling visualization

* [ ] time estimates

* [ ] color smoothing optional

* [ ] publish `0.1.0`

* [ ] wasm based editor

* [ ] drawing algorithms: statistical methods (mean, avg, etc.)

* [ ] `julia set`: field lines

* [ ] progressive iteration visualization as GIF

* [ ] newton fractals

* [ ] fractal flame


### MAYBEDO

* [ ] support for alpha channel

* [ ] greyscale images

* [ ] more gradients

* [ ] color map: stops for colors

* [ ] arbitrary precision types for higher zoom levels
