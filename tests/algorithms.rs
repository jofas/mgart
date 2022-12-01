use mgart::buddhabrot::Buddhabrot;
use mgart::julia_set::JuliaSet;
use mgart::util::coloring::ColorMap1d;
use mgart::util::post_processing::PostProcessing;
use mgart::util::sampler::Sampler;
use mgart::util::ComplexNumber;

#[test]
fn buddhabrot() {
  let buddhabrot = Buddhabrot::new(
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
  );

  drop(buddhabrot.creator().create());
}

#[test]
fn julia_set() {
  let julia_set = JuliaSet::new(
    1920,
    1024,
    ComplexNumber::Cartesian { re: 0., im: 0. },
    1.,
    20,
    None,
    ColorMap1d::default(),
    None,
  );

  drop(julia_set.creator().create());
}
