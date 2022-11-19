use serde::{Deserialize, Serialize};

use num_complex::Complex64;

pub mod coloring;
pub mod errors;
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

pub fn print_progress(i: u64, n: u64, interval: u64) {
  if i % interval == interval - 1 || i == n - 1 {
    let p = i as f64 / n as f64 * 100.;
    print!("{}/{} iterations done ({:.2}%)\r", i + 1, n, p);
  }
}
