use num_complex::Complex64;

use rand::random;

use std::f64::consts::PI;
use std::marker::PhantomData;

pub trait Sampler {
  type Output;

  fn sample(&self) -> Self::Output;
}

pub struct Uniform<T> {
  phantom: PhantomData<T>,
}

impl<T> Uniform<T> {
  pub fn new() -> Self {
    Self {
      phantom: PhantomData,
    }
  }
}

impl Sampler for Uniform<Complex64> {
  type Output = Complex64;

  fn sample(&self) -> Self::Output {
    random_complex()
  }
}

/// Kernel Density Estimation.
///
pub struct KDE<T, K> {
  elems: Vec<T>,
  probabilities: Vec<f64>,
  kernel: K,
}

impl<T, K: Fn(&T) -> T> KDE<T, K> {
  pub fn new(samples: Vec<(T, f64)>, kernel: K) -> Self {
    let (elems, probabilities) = samples.into_iter().unzip();

    Self {
      elems: elems,
      probabilities: probabilities,
      kernel,
    }
  }

  pub fn average_probability(&self) -> f64 {
    self.probabilities.iter().sum::<f64>() / self.elems.len() as f64
  }
}

impl<T, K: Fn(&T) -> T> Sampler for KDE<T, K> {
  type Output = T;

  fn sample(&self) -> Self::Output {
    loop {
      let idx = (random::<f64>() * self.elems.len() as f64) as usize;

      if random::<f64>() <= self.probabilities[idx] {
        return (self.kernel)(&self.elems[idx]);
      }
    }
  }
}

/// Random complex number with a radius between `[0, 2]` from the
/// origin.
///
pub fn random_complex() -> Complex64 {
  Complex64::from_polar(
    2. * random::<f64>(),
    2. * PI * random::<f64>(),
  )
}
