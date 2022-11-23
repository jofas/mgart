use serde::{Deserialize, Serialize};

use anyhow::Result;

use crate::util::gradient::Gradient;

pub mod clahe;
pub mod smoothing;

use clahe::CLAHE;
use smoothing::Smoothing;

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "process")]
pub enum PostProcessing {
  Normalize,
  Clamp { min: f64, max: f64 },
  ClampAndNormalize { min: f64, max: f64 },
  Gradient(Gradient),
  Smoothing(Smoothing),
  Clahe(CLAHE),
}

impl PostProcessing {
  /// Applies the [`PostProcessing`] algorithm to `buffer`.
  ///
  /// # Errors
  ///
  /// An algorithm may fail, if it's configuration is faulty.
  /// Otherwise, an error running a [`PostProcessing`] algorithm is
  /// a bug.
  ///
  pub fn apply(
    &self,
    buffer: &mut [f64],
    width: usize,
    height: usize,
  ) -> Result<()> {
    match self {
      Self::Normalize => {
        let (min, max) = min_max(buffer);

        for v in buffer {
          *v = (*v - min) / (max - min);
        }
      }
      Self::Clamp { min, max } => {
        for v in buffer {
          *v = v.clamp(*min, *max);
        }
      }
      Self::ClampAndNormalize { min, max } => {
        for v in buffer {
          *v = (v.clamp(*min, *max) - min) / (max - min);
        }
      }
      Self::Gradient(g) => {
        for v in buffer {
          *v = g.apply(*v);
        }
      }
      Self::Smoothing(s) => {
        s.smooth(buffer, width, height);
      }
      Self::Clahe(c) => {
        c.apply(buffer, width, height)?;
      }
    }

    Ok(())
  }
}

fn min_max(v: &[f64]) -> (f64, f64) {
  let mut max = &0.;
  let mut min = &f64::MAX;

  for x in v {
    if x < min {
      min = x;
    }

    if x > max {
      max = x;
    }
  }

  (*min, *max)
}
