use criterion::{criterion_group, criterion_main, Criterion};

use std::time::Duration;

use mgart::buddhabrot::Buddhabrot;
use mgart::util::coloring::ColorMap1d;
use mgart::util::sampler::Sampler;
use mgart::util::ComplexNumber;

pub fn buddhabrot_kde_sampler(c: &mut Criterion) {
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
            population: 10_000_000,
            p_min: 0.01,
            pre_sampler: Box::new(Sampler::UniformPolar { r: 3. }),
        },
        vec![],
    );

    let creator = buddhabrot.creator();

    c.bench_function("buddhabrot kde sampler", |b| {
        b.iter(|| creator.clone().sampler());
    });
}

criterion_group!(
  name = benches;
  config = Criterion::default()
    .sample_size(10)
    .measurement_time(Duration::new(10, 0));
  targets = buddhabrot_kde_sampler
);

criterion_main!(benches);
