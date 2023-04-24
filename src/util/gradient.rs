use serde::{Deserialize, Serialize};

use std::f64::consts::PI;

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "type")]
pub enum Gradient {
  Linear { factor: f64 },
  Sin { factor: f64 },
  Inverted { gradient: Box<Gradient> },
  Wave { factor: f64 },
  Exp { exponent: f64 },
  SinExp { factor: f64 },
  Log { factor: f64 },
  Tanh { factor: f64 },
  SinRamp { factor: f64, amplitude: f64 },
  Discrete { gradient: Box<Gradient>, n: u32 },
  Smoothstep { order: i32 },
  // TODO: b-spline
}

impl Gradient {
  #[must_use]
  pub fn apply(&self, f: f64) -> f64 {
    match self {
      Self::Linear { factor } => {
        let f = f.clamp(0., 1.);

        let res = (f * factor).fract();

        if f <= f64::EPSILON {
          0.
        } else if res <= f64::EPSILON {
          1.
        } else {
          (f * factor).fract()
        }
      }
      Self::Sin { factor } => (f * factor * PI).sin() / 2. + 0.5,
      Self::Inverted { gradient } => 1. - gradient.apply(f),
      Self::Wave { factor } => {
        let f = (f * factor).fract();

        if f <= 0.5 {
          f
        } else {
          1. - f
        }
      }
      Self::Exp { exponent } => f.powf(*exponent),
      Self::SinExp { factor } => (f * factor * PI).exp().sin().abs(),
      Self::Log { factor } => {
        (f * factor + 1.).ln() / (factor + 1.).ln()
      }
      Self::Tanh { factor } => (f * factor).tanh() / factor.tanh(),
      Self::SinRamp { factor, amplitude } => {
        (f + amplitude * (f * factor * PI)).clamp(0., 1.)
      }
      Self::Discrete { gradient, n } => {
        let n = f64::from(*n);
        let f = gradient.apply(f).clamp(0., 1.);

        let res = (f * n).floor() / (n - 1.);

        if (1. - f).abs() <= f64::EPSILON {
          res - 1. / (n - 1.).max(1.)
        } else {
          res
        }
      }
      Self::Smoothstep { order } => {
        let f = f.clamp(0., 1.);
        (0..=*order).fold(0., |acc, n| {
          let p1 = f64::from(Self::pascal_triangle(-order - 1, n));

          let p2 = f64::from(Self::pascal_triangle(
            2 * order + 1,
            order - n,
          ));

          acc + p1 * p2 * f.powi(order + n + 1)
        })
      }
    }
  }

  fn pascal_triangle(a: i32, b: i32) -> i32 {
    (0..b).fold(1, |acc, i| acc * (a - i) / (i + 1))
  }
}

impl Default for Gradient {
  fn default() -> Self {
    Self::Linear { factor: 1. }
  }
}

#[cfg(test)]
mod tests {
  use super::Gradient;

  #[test]
  fn linear_gradient() {
    let g = Gradient::Linear { factor: 1. };

    assert!((g.apply(0.) - 0.).abs() <= f64::EPSILON);
    assert!((g.apply(0.5) - 0.5).abs() <= f64::EPSILON);
    assert!((g.apply(1.) - 1.).abs() <= f64::EPSILON);

    let g = Gradient::Linear { factor: 2. };

    assert!((g.apply(0.) - 0.).abs() <= f64::EPSILON);
    assert!((g.apply(0.25) - 0.5).abs() <= f64::EPSILON);
    assert!((g.apply(0.5) - 1.).abs() <= f64::EPSILON);
    assert!((g.apply(0.75) - 0.5).abs() <= f64::EPSILON);
    assert!((g.apply(1.) - 1.).abs() <= f64::EPSILON);
  }

  #[test]
  fn discrete_gradient() {
    let g = Gradient::Discrete {
      n: 2,
      gradient: Box::new(Gradient::Linear { factor: 1. }),
    };

    assert!((g.apply(0.) - 0.).abs() <= f64::EPSILON);
    assert!((g.apply(0.25) - 0.).abs() <= f64::EPSILON);
    assert!((g.apply(0.5) - 1.).abs() <= f64::EPSILON);
    assert!((g.apply(0.75) - 1.).abs() <= f64::EPSILON);
    assert!((g.apply(1.) - 1.).abs() <= f64::EPSILON);

    let g = Gradient::Discrete {
      n: 3,
      gradient: Box::new(Gradient::Linear { factor: 1. }),
    };

    assert!((g.apply(0.) - 0.).abs() <= f64::EPSILON);
    assert!((g.apply(0.33) - 0.).abs() <= f64::EPSILON);
    assert!((g.apply(0.34) - 0.5).abs() <= f64::EPSILON);
    assert!((g.apply(0.66) - 0.5).abs() <= f64::EPSILON);
    assert!((g.apply(0.67) - 1.).abs() <= f64::EPSILON);
    assert!((g.apply(1.) - 1.).abs() <= f64::EPSILON);

    let g = Gradient::Discrete {
      n: 4,
      gradient: Box::new(Gradient::Linear { factor: 1. }),
    };

    assert!((g.apply(0.) - 0.).abs() <= f64::EPSILON);
    assert!((g.apply(0.24) - 0.).abs() <= f64::EPSILON);
    assert!(
      ((g.apply(0.25) * 100.).floor() - 33.).abs() <= f64::EPSILON
    );
    assert!(
      ((g.apply(0.49) * 100.).floor() - 33.).abs() <= f64::EPSILON
    );
    assert!(
      ((g.apply(0.5) * 100.).floor() - 66.).abs() <= f64::EPSILON
    );
    assert!(
      ((g.apply(0.74) * 100.).floor() - 66.).abs() <= f64::EPSILON
    );
    assert!((g.apply(0.75) - 1.).abs() <= f64::EPSILON);
    assert!((g.apply(1.) - 1.).abs() <= f64::EPSILON);
  }
}
