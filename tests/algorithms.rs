use mgart::buddhabrot::Buddhabrot;
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
    200,
    None,
    ColorMap1d::default(),
    2.,
    100_000,
    Sampler::UniformPolar { r: 2. },
    vec![PostProcessing::Normalize],
  );

  buddhabrot.creator().create();
}
