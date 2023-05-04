use criterion::{criterion_group, criterion_main, Criterion};

use rand::random;

use map_macro::vec_no_clone;

use std::time::Duration;

use mgart::util::frame::Frame;
use mgart::util::post_processing::clahe::CLAHE;
use mgart::util::post_processing::smoothing::non_local_means::NonLocalMeans;

pub fn clahe(c: &mut Criterion) {
    let clahe = CLAHE::new(60, 256, 1100, 1100);

    let width = 11_000;
    let height = 11_000;

    let buffer = vec_no_clone![random::<f64>(); width * height];

    let mut frame = Frame::new(buffer, width, height);

    c.bench_function("clahe", |b| {
        b.iter(|| clahe.apply(&mut frame));
    });
}

pub fn nlm(c: &mut Criterion) {
    let nlm = NonLocalMeans::new(7, 21, 1e-4);

    let width = 700;
    let height = 700;

    let buffer = vec_no_clone![random::<f64>(); width * height];

    let mut frame = Frame::new(buffer, width, height);

    c.bench_function("nlm", |b| {
        b.iter(|| nlm.smooth(&mut frame));
    });
}

criterion_group!(
  name = benches;
  config = Criterion::default()
    .sample_size(10)
    .measurement_time(Duration::new(10, 0));
  targets = clahe, nlm
);

criterion_main!(benches);
