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

  /// Generates a uniformly random complex number with a radius between
  /// `[0, 2]` from the origin.
  ///
  fn sample(&self) -> Self::Output {
    Complex64::from_polar(
      2. * random::<f64>(),
      2. * PI * random::<f64>(),
    )
  }
}

/// Kernel Density Estimation.
///
pub struct KDE<T, K> {
  elems: Vec<T>,
  kernel: K,
}

impl<T, K: Fn(&T) -> T> KDE<T, K> {
  pub fn new(samples: Vec<T>, kernel: K) -> Self {
    Self {
      elems: samples,
      kernel,
    }
  }
}

impl<T, K: Fn(&T) -> T> Sampler for KDE<T, K> {
  type Output = T;

  fn sample(&self) -> Self::Output {
    let idx = (random::<f64>() * self.elems.len() as f64) as usize;
    return (self.kernel)(&self.elems[idx]);
  }
}

/// Weighted Kernel Density Estimation.
///
pub struct WeightedKDE<T, K> {
  elems: Vec<T>,
  probabilities: Vec<f64>,
  kernel: K,
}

impl<T, K: Fn(&T) -> T> WeightedKDE<T, K> {
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

impl<T, K: Fn(&T) -> T> Sampler for WeightedKDE<T, K> {
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
