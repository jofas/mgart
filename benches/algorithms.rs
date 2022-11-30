use criterion::{criterion_group, criterion_main, Criterion};

use std::time::Duration;

use mgart::buddhabrot::Buddhabrot;
use mgart::util::coloring::ColorMap1d;
use mgart::util::sampler::Sampler;
use mgart::util::ComplexNumber;

pub fn buddhabrot(c: &mut Criterion) {
  let buddhabrot = Buddhabrot::new(
    1920,
    1024,
    ComplexNumber::Cartesian { re: 0., im: 0. },
    1.,
    200,
    None,
    ColorMap1d::default(),
    2.,
    1_000_000,
    Sampler::UniformPolar { r: 2. },
    vec![],
  );

  let creator = buddhabrot.creator();

  c.bench_function("buddhabrot", |b| {
    b.iter(|| creator.clone().create());
  });
}

criterion_group!(
  name = benches;
  config = Criterion::default()
    .sample_size(10)
    .measurement_time(Duration::new(10, 0));
  targets = buddhabrot
);

criterion_main!(benches);
