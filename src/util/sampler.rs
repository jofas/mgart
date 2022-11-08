use serde::{Deserialize, Serialize};

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use num_complex::Complex64;

use rand::random;

use std::f64::consts::PI;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::util::print_progress;

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "type")]
pub enum Sampler {
  /// Uniform random distribution generating values in the range of
  /// `[0,h]`.
  ///
  Uniform { h: f64 },
  /// Uniform random distribution generating values with a radius
  /// in `[0,r]` from the origin of a space, e.g. `(0,0)` in an
  /// unspecified 2d space or `0 + 0i` in the complex space.
  ///
  UniformPolar { r: f64 },
  /// Kernel Density Estimation.
  ///
  /// *NOTE:* `kernel` and `pre_sampler` must eventually resolve to a
  /// sampler that does not require pre-sampling, e.g.
  /// [Self::Uniform] or [Self::UniformPolar].
  ///
  KernelDensityEstimation {
    kernel: Box<Sampler>,
    population: u64,
    p_min: f64,
    pre_sampler: Box<Sampler>,
  },
  WeightedKernelDensityEstimation {
    kernel: Box<Sampler>,
    population: u64,
    p_min: f64,
    bin_count: u64,
    pre_sampler: Box<Sampler>,
  },
}

impl Sampler {
  /// Creates a [Distribution] from this instance of [Sampler].
  ///
  /// The `true_probability` function is needed for the
  /// [Self::KernelDensityEstimation] to create samples from the
  /// true probability function.
  ///
  pub fn distribution<
    T: Sync + Send,
    P: Fn(&T) -> f64 + Sync + Send,
  >(
    self,
    true_probability: &P,
  ) -> Distribution<T>
  where
    Distribution<T>: Sampling<Space = T>,
  {
    match self {
      Self::Uniform { h } => {
        Distribution::Uniform(Uniform::<T>::new(h))
      }
      Self::UniformPolar { r } => {
        Distribution::UniformPolar(UniformPolar::<T>::new(r))
      }
      Self::KernelDensityEstimation {
        kernel,
        population,
        p_min,
        pre_sampler,
      } => {
        let (elems, _) = Self::pre_sample(
          pre_sampler.distribution(true_probability),
          true_probability,
          population,
          p_min,
        );

        Distribution::KernelDensityEstimation(KDE::<T>::new(
          elems,
          Box::new(kernel.distribution(true_probability)),
        ))
      }
      Self::WeightedKernelDensityEstimation {
        kernel,
        population,
        p_min,
        bin_count,
        pre_sampler,
      } => {
        let (elems, probabilities) = Self::pre_sample(
          pre_sampler.distribution(true_probability),
          true_probability,
          population,
          p_min,
        );

        Distribution::WeightedKernelDensityEstimation(
          WeightedKDE::<T>::new(
            elems,
            probabilities,
            bin_count,
            Box::new(kernel.distribution(true_probability)),
          ),
        )
      }
    }
  }

  fn pre_sample<T: Sync + Send, P: Fn(&T) -> f64 + Sync + Send>(
    pre_sampler: Distribution<T>,
    true_probability: &P,
    population: u64,
    p_min: f64,
  ) -> (Vec<T>, Vec<f64>)
  where
    Distribution<T>: Sampling<Space = T>,
    Vec<T>: Extend<T>,
  {
    println!("initializing kde population");

    let processed_samples = AtomicU64::new(0);

    let res = (0..population)
      .into_par_iter()
      .filter_map(|_| {
        let sample = pre_sampler.sample();

        let p = true_probability(&sample);

        let ps = processed_samples.fetch_add(1, Ordering::SeqCst);
        print_progress(ps, population, 2500);

        if p >= p_min {
          Some((sample, p))
        } else {
          None
        }
      })
      .unzip();

    println!("\ninitializing kde population done");

    res
  }
}

pub trait Sampling {
  type Space;

  fn sample(&self) -> Self::Space;
  fn kernel_sample(&self, center: &Self::Space) -> Self::Space;
}

pub enum Distribution<T> {
  Uniform(Uniform<T>),
  UniformPolar(UniformPolar<T>),
  KernelDensityEstimation(KDE<T>),
  WeightedKernelDensityEstimation(WeightedKDE<T>),
}

impl Sampling for Distribution<Complex64> {
  type Space = Complex64;

  fn sample(&self) -> Self::Space {
    match self {
      Self::Uniform(u) => u.sample(),
      Self::UniformPolar(up) => up.sample(),
      Self::KernelDensityEstimation(kde) => kde.sample(),
      Self::WeightedKernelDensityEstimation(wkde) => wkde.sample(),
    }
  }

  fn kernel_sample(&self, center: &Self::Space) -> Self::Space {
    match self {
      Self::Uniform(u) => u.kernel_sample(center),
      Self::UniformPolar(up) => up.kernel_sample(center),
      Self::KernelDensityEstimation(kde) => kde.kernel_sample(center),
      Self::WeightedKernelDensityEstimation(wkde) => {
        wkde.kernel_sample(center)
      }
    }
  }
}

pub struct Uniform<T> {
  h: f64,
  phantom: PhantomData<T>,
}

impl<T> Uniform<T> {
  pub fn new(h: f64) -> Self {
    Self {
      h,
      phantom: PhantomData,
    }
  }
}

impl Sampling for Uniform<Complex64> {
  type Space = Complex64;

  /// Generates a uniformly random complex number with real and
  /// imaginary parts ranging from `[0, h]`.
  ///
  fn sample(&self) -> Self::Space {
    Complex64::new(random::<f64>() * self.h, random::<f64>() * self.h)
  }

  /// Generates a uniformly random complex number with real and
  /// imaginary parts ranging from `[center - h/2, center + h/2]`.
  ///
  fn kernel_sample(&self, center: &Self::Space) -> Self::Space {
    let re = (random::<f64>() - 0.5) * self.h;
    let im = (random::<f64>() - 0.5) * self.h;

    Complex64::new(center.re + re, center.im + im)
  }
}

pub struct UniformPolar<T> {
  r: f64,
  phantom: PhantomData<T>,
}

impl<T> UniformPolar<T> {
  pub fn new(r: f64) -> Self {
    Self {
      r,
      phantom: PhantomData,
    }
  }
}

impl Sampling for UniformPolar<Complex64> {
  type Space = Complex64;

  /// Generates a uniformly random complex number with a radius `r`
  /// from the origin of the complex space `0 + 0i`.
  ///
  fn sample(&self) -> Self::Space {
    Complex64::from_polar(
      self.r * random::<f64>(),
      2. * PI * random::<f64>(),
    )
  }

  fn kernel_sample(&self, center: &Self::Space) -> Self::Space {
    center + self.sample()
  }
}

/// Kernel Density Estimation.
///
pub struct KDE<T> {
  elems: Vec<T>,
  kernel: Box<Distribution<T>>,
}

impl<T> KDE<T> {
  pub fn new(elems: Vec<T>, kernel: Box<Distribution<T>>) -> Self {
    Self { elems, kernel }
  }
}

impl Sampling for KDE<Complex64> {
  type Space = Complex64;

  fn sample(&self) -> Self::Space {
    let idx = (random::<f64>() * self.elems.len() as f64) as usize;
    return self.kernel.kernel_sample(&self.elems[idx]);
  }

  fn kernel_sample(&self, center: &Self::Space) -> Self::Space {
    center + self.sample()
  }
}

/// Weighted Kernel Density Estimation.
///
pub struct WeightedKDE<T> {
  elems: Vec<T>,
  probabilities: Vec<f64>,
  kernel: Box<Distribution<T>>,
}

impl<T> WeightedKDE<T> {
  pub fn new(
    elems: Vec<T>,
    probabilities: Vec<f64>,
    bin_count: u64,
    kernel: Box<Distribution<T>>,
  ) -> Self {
    // TODO: create bins
    //       for each bin derive its probability
    //       create list where we sample from each bin according to
    //       its probability
    //       then sample from weighted kde like from kde
    Self {
      elems,
      probabilities,
      kernel,
    }
  }
}

impl Sampling for WeightedKDE<Complex64> {
  type Space = Complex64;

  fn sample(&self) -> Self::Space {
    loop {
      let idx = (random::<f64>() * self.elems.len() as f64) as usize;

      if random::<f64>() <= self.probabilities[idx] {
        return self.kernel.kernel_sample(&self.elems[idx]);
      }
    }
  }

  fn kernel_sample(&self, center: &Self::Space) -> Self::Space {
    center + self.sample()
  }
}
