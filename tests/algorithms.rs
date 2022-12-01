use std::fs::remove_file;

use mgart::buddhabrot::Buddhabrot;
use mgart::debug::ColorMap1dRenderer;
use mgart::julia_set::JuliaSet;
use mgart::util::coloring::ColorMap1d;
use mgart::util::post_processing::PostProcessing;
use mgart::util::sampler::Sampler;
use mgart::util::ComplexNumber;
use mgart::{Algorithm, Algorithms};

#[test]
fn buddhabrot() {
  let filename = "test_buddhabrot.png";

  let buddhabrot = Algorithm::buddhabrot(
    Buddhabrot::new(
      1920,
      1024,
      ComplexNumber::Cartesian { re: 0., im: 0. },
      1.,
      20,
      None,
      ColorMap1d::default(),
      2.,
      1000,
      Sampler::KernelDensityEstimation {
        weighted: true,
        kernel: Box::new(Sampler::Uniform { h: 2e-2 }),
        population: 1000,
        p_min: 0.01,
        pre_sampler: Box::new(Sampler::UniformPolar { r: 3. }),
      },
      vec![PostProcessing::Normalize],
    ),
    filename.to_owned(),
  );

  buddhabrot.create();

  remove_file(filename).unwrap();
}

#[test]
fn julia_set_mandelbrot() {
  let filename = "test_julia_set_mandelbrot.png";

  let julia_set = Algorithm::julia_set(
    JuliaSet::new(
      256,
      256,
      ComplexNumber::Cartesian { re: 0., im: 0. },
      1.,
      20,
      None,
      ColorMap1d::default(),
      None,
    ),
    filename.to_owned(),
  );

  julia_set.create();

  remove_file(filename).unwrap();
}

#[test]
fn julia_set_julia_set() {
  env_logger::Builder::from_env(
    env_logger::Env::default().default_filter_or("info"),
  )
  .init();

  let filename = "test_julia_set_julia_set.png";

  let julia_set = Algorithm::julia_set(
    JuliaSet::new(
      256,
      256,
      ComplexNumber::Cartesian { re: 0., im: 0. },
      1.,
      20,
      None,
      ColorMap1d::default(),
      Some(ComplexNumber::Polar { r: 0., theta: 0. }),
    ),
    filename.to_owned(),
  );

  julia_set.create();

  remove_file(filename).unwrap();
}

#[test]
fn color_map_1d_renderer() {
  let filename = "test_color_map_1d_renderer.png";

  let cmr = Algorithms::new(vec![Algorithm::color_map_1d_renderer(
    ColorMap1dRenderer::new(256, 256, ColorMap1d::default()),
    filename.to_owned(),
  )]);

  cmr.create();

  remove_file(filename).unwrap();
}
