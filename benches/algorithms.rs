use criterion::{criterion_group, criterion_main, Criterion};

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
    20,
    None,
    ColorMap1d::default(),
    2.,
    1000,
    Sampler::UniformPolar { r: 2. },
    vec![],
  );

  c.bench_function("buddhabrot", |b| {
    b.iter(|| buddhabrot.clone().create());
  });
}

criterion_group!(benches, buddhabrot);
criterion_main!(benches);
