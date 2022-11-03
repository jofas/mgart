use serde::{Deserialize, Serialize};

use num_complex::Complex64;

pub mod coloring;
pub mod gradient;
pub mod post_processing;
pub mod sampler;
pub mod viewport;

/// Representation of a complex number.
///
/// This is intended to be used as means for parsing user input,
/// not for doing calculations.
/// So [ComplexNumber] does not implement any math operations,
/// but supports the conversion to [Complex64].
///
#[derive(Serialize, Deserialize, Clone, Copy)]
#[serde(untagged)]
pub enum ComplexNumber {
  Cartesian { re: f64, im: f64 },
  Polar { r: f64, theta: f64 },
}

impl Into<Complex64> for ComplexNumber {
  fn into(self) -> Complex64 {
    match self {
      ComplexNumber::Cartesian { re, im } => Complex64::new(re, im),
      ComplexNumber::Polar { r, theta } => {
        Complex64::from_polar(r, theta)
      }
    }
  }
}

impl Into<Complex64> for &ComplexNumber {
  fn into(self) -> Complex64 {
    (*self).into()
  }
}

pub fn print_progress(i: u64, n: u64, interval: u64) {
  if i % interval == interval - 1 || i == n - 1 {
    let p = i as f64 / n as f64 * 100.;
    print!("{}/{} iterations done ({:.2}%)\r", i + 1, n, p);
  }
}
