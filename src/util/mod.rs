use serde::{Deserialize, Serialize};

use num::complex::Complex64;

use log::{log_enabled, Level};

use num::cast;

use std::sync::atomic::{AtomicU64, Ordering};

pub mod coloring;
pub mod frame;
pub mod gradient;
pub mod post_processing;
pub mod sampler;
pub mod viewport;

/// Representation of a complex number.
///
/// This is intended to be used as means for parsing user input,
/// not for doing calculations.
/// So [`ComplexNumber`] does not implement any math operations,
/// but supports the conversion to [Complex64].
///
#[derive(Serialize, Deserialize, Clone, Copy)]
#[serde(untagged)]
pub enum ComplexNumber {
  Cartesian { re: f64, im: f64 },
  Polar { r: f64, theta: f64 },
}

impl From<&ComplexNumber> for Complex64 {
  fn from(cn: &ComplexNumber) -> Self {
    match cn {
      ComplexNumber::Cartesian { re, im } => Self::new(*re, *im),
      ComplexNumber::Polar { r, theta } => {
        Self::from_polar(*r, *theta)
      }
    }
  }
}

impl From<ComplexNumber> for Complex64 {
  fn from(cn: ComplexNumber) -> Self {
    Self::from(&cn)
  }
}

// TODO: dz is derivative of z in the iteration sequence ... if I were
//       to look for finite attractors for different functions, I'd
//       need to generalize this
#[must_use]
pub fn finite_attractor(
  z0: Complex64,
  c: Complex64,
  p: usize,
) -> Option<(Complex64, Complex64)> {
  let mut zz = z0;

  for _ in 0..64 {
    let mut z = zz;
    let mut dz = Complex64::new(1., 0.);

    for _ in 0..p {
      dz = 2. * z * dz;
      z = z.powi(2) + c;
    }

    let zz_new = zz - (z - zz) / (dz - 1.);

    if (zz_new - zz).norm_sqr() <= 1e-20 {
      return Some((z, dz));
    }

    zz = zz_new;
  }

  None
}

/// Thread-safe progress printer for printing the number of iterations
/// done by a computationally expensive loop.
///
pub enum ProgressPrinter {
  Info(InfoProgressPrinter),
  Disabled,
}

impl ProgressPrinter {
  #[must_use]
  pub fn new(n: u64, interval: u64) -> Self {
    if log_enabled!(Level::Info) {
      Self::Info(InfoProgressPrinter::new(n, interval))
    } else {
      Self::Disabled
    }
  }

  pub fn increment(&self) {
    if let Self::Info(info) = self {
      info.increment();
    }
  }
}

pub struct InfoProgressPrinter {
  counter: AtomicU64,
  n: u64,
  n_f64: f64,
  interval: u64,
}

impl InfoProgressPrinter {
  /// Creates a new instance of [`InfoProgressPrinter`].
  ///
  /// # Panics
  ///
  /// Panics, if [`n`] overflows the 53 bit mantissa of [f64].
  ///
  #[must_use]
  pub fn new(n: u64, interval: u64) -> Self {
    Self {
      counter: AtomicU64::new(0),
      n,
      n_f64: cast::<_, f64>(n).unwrap(),
      interval,
    }
  }

  /// Increments the counter, printing a status update if the counter
  /// equals a multiple of [`interval`] or [`n`].
  ///
  /// # Panics
  ///
  /// Panics, if the counter overflows the 53 bit mantissa of [f64].
  ///
  pub fn increment(&self) {
    let i = self.counter.fetch_add(1, Ordering::SeqCst);

    if i % self.interval == self.interval - 1 || i == self.n - 1 {
      let p = cast::<_, f64>(i).unwrap() / self.n_f64 * 100.;
      print!("{}/{} iterations done ({:.2}%)\r", i + 1, self.n, p);
    }
  }
}
