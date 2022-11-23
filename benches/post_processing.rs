use criterion::{criterion_group, criterion_main, Criterion};

use rand::random;

use map_macro::vec_no_clone;

use mgart::util::post_processing::clahe::CLAHE;

pub fn clahe(c: &mut Criterion) {
  let clahe = CLAHE::new(60, 256, 64, 64);

  let width = 640;
  let height = 640;

  let mut buffer = vec_no_clone![random::<f64>(); width * height];

  c.bench_function("clahe", |b| {
    b.iter(|| clahe.apply(&mut buffer, width, height));
  });
}

criterion_group!(benches, clahe);
criterion_main!(benches);
