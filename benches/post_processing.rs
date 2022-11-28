use criterion::{criterion_group, criterion_main, Criterion};

use rand::random;

use map_macro::vec_no_clone;

use mgart::util::frame::Frame;
use mgart::util::post_processing::clahe::CLAHE;
use mgart::util::post_processing::smoothing::non_local_means::NonLocalMeans;

pub fn clahe(c: &mut Criterion) {
  let clahe = CLAHE::new(60, 256, 64, 64);

  let buffer = vec_no_clone![random::<f64>(); 1024 * 1024];

  let mut frame = Frame::new(buffer, 1024, 1024);

  c.bench_function("clahe", |b| {
    b.iter(|| clahe.apply(&mut frame));
  });
}

pub fn nlm(c: &mut Criterion) {
  let nlm = NonLocalMeans::new(7, 21, 1e-4);

  let buffer = vec_no_clone![random::<f64>(); 256 * 256];

  let mut frame = Frame::new(buffer, 256, 256);

  c.bench_function("nlm", |b| {
    b.iter(|| nlm.smooth(&mut frame));
  });
}

criterion_group!(benches, clahe, nlm);
criterion_main!(benches);
