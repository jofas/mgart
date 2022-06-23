use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

use num_complex::Complex64;

use map_macro::vec_no_clone;

use rand::random;
use rand_distr::{Distribution, Normal};

use std::f64::consts::PI;

pub mod colors;

use colors::{Color, LCH, RGB};

/// Representation of a complex number.
///
/// This is intended to be used as means for parsing user input,
/// not for doing calculations.
/// So [ComplexNumber] does not implement any math operations,
/// but supports the conversion to [Complex64].
///
#[derive(Serialize, Deserialize)]
#[serde(untagged)]
pub enum ComplexNumber {
  Cartesian { re: f64, im: f64 },
  Polar { r: f64, theta: f64 },
}

impl Into<Complex64> for &ComplexNumber {
  fn into(self) -> Complex64 {
    match self {
      ComplexNumber::Cartesian { re, im } => Complex64::new(*re, *im),
      ComplexNumber::Polar { r, theta } => {
        Complex64::from_polar(*r, *theta)
      }
    }
  }
}

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
  pub fn apply_to(&self, f: f64) -> f64 {
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
      Self::Inverted { gradient } => 1. - gradient.apply_to(f),
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
        let n = *n as f64;
        let f = gradient.apply_to(f).clamp(0., 1.);

        let res = (f * n).floor() / (n - 1.);

        if (1. - f).abs() <= f64::EPSILON {
          res - 1. / (n - 1.).max(1.)
        } else {
          res
        }
      }
      Self::Smoothstep { order } => {
        let f = f.clamp(0., 1.);
        (0..=*order).into_iter().fold(0., |acc, n| {
          acc
            + Self::pascal_triangle(-order - 1, n) as f64
              * Self::pascal_triangle(2 * order + 1, order - n) as f64
              * f.powi(order + n + 1)
        })
      }
    }
  }

  fn pascal_triangle(a: i32, b: i32) -> i32 {
    (0..b).into_iter().fold(1, |acc, i| acc * (a - i) / (i + 1))
  }
}

#[derive(Serialize, Deserialize, DisplayAsJson)]
#[serde(from = "ColorMap1dDeserializer")]
pub struct ColorMap1d {
  map: Vec<LCH>,
  gradient: Gradient,
}

impl ColorMap1d {
  pub fn new(map: Vec<Color>, gradient: Gradient) -> Self {
    let map = if map.len() >= 2 {
      map
    } else if map.len() == 1 {
      vec![Color::WHITE, map[0]]
    } else {
      vec![Color::WHITE, Color::BLACK]
    };

    let map: Vec<LCH> = map.into_iter().map(|c| c.lch()).collect();

    Self { map, gradient }
  }

  pub fn color(&self, f: f64) -> RGB {
    let f = self.gradient.apply_to(f);

    if (1.0 - f).abs() <= f64::EPSILON {
      return self.map[self.map.len() - 1].rgb();
    }

    let interval = f * (self.map.len() - 1) as f64;
    let pos = interval.fract() as f64;

    let c1 = &self.map[interval as usize];
    let c2 = &self.map[interval as usize + 1];

    c1.interpolate(&c2, pos).rgb()
  }
}

impl Default for ColorMap1d {
  fn default() -> Self {
    Self::new(vec![], Gradient::Linear { factor: 1. })
  }
}

#[derive(Deserialize)]
struct ColorMap1dDeserializer {
  map: Vec<Color>,
  gradient: Gradient,
}

impl From<ColorMap1dDeserializer> for ColorMap1d {
  fn from(cm: ColorMap1dDeserializer) -> Self {
    Self::new(cm.map, cm.gradient)
  }
}

pub struct GaussianKDE {
  len: usize,
  elems: Vec<Complex64>,
  probabilities: Vec<f64>,
  distribution: Normal<f64>,
}

impl GaussianKDE {
  pub fn new(len: usize, h: f64) -> Self {
    Self {
      len,
      elems: vec_no_clone![random_complex(); len],
      probabilities: vec![0.; len],
      distribution: Normal::new(0., h).unwrap(),
    }
  }

  pub fn sample(&self) -> Complex64 {
    let idx = (random::<f64>() * self.elems.len() as f64) as usize;

    let re = self.distribution.sample(&mut rand::thread_rng());
    let im = self.distribution.sample(&mut rand::thread_rng());

    Complex64::new(self.elems[idx].re + re, self.elems[idx].im + im)
  }

  pub fn update(&mut self, c: Complex64, p: f64) {
    // TODO: don't allow similar points
    //
    // ... how?
    //
    // get nn, if better than nn, replace, else, keep nn
    //

    let idx = (random::<f64>() * self.len as f64) as usize;

    if self.probabilities[idx] < p {
      self.elems[idx] = c;
      self.probabilities[idx] = p;
    }
  }

  pub fn average_probability(&self) -> f64 {
    self.probabilities.iter().sum::<f64>() / self.len as f64
  }
}

pub fn in_viewport(
  x: f64,
  y: f64,
  x_min: f64,
  x_max: f64,
  y_min: f64,
  y_max: f64,
) -> bool {
  x_min <= x && x < x_max && y_min <= y && y < y_max
}

pub fn grid_pos(
  x: f64,
  y: f64,
  x_min: f64,
  x_max: f64,
  y_min: f64,
  y_max: f64,
  delta_x: f64,
  delta_y: f64,
) -> Option<(usize, usize)> {
  if in_viewport(x, y, x_min, x_max, y_min, y_max) {
    let x = ((x - x_min) / delta_x) as usize;
    let y = ((y - y_min) / delta_y) as usize;

    Some((x, y))
  } else {
    None
  }
}

pub fn random_complex() -> Complex64 {
  Complex64::from_polar(
    2. * random::<f64>(),
    2. * PI * random::<f64>(),
  )
}

#[cfg(test)]
mod tests {
  use super::{grid_pos, Gradient};

  #[test]
  fn linear_gradient() {
    let g = Gradient::Linear { factor: 1. };

    assert_eq!(g.apply_to(0.), 0.);
    assert_eq!(g.apply_to(0.5), 0.5);
    assert_eq!(g.apply_to(1.), 1.);

    let g = Gradient::Linear { factor: 2. };

    assert_eq!(g.apply_to(0.), 0.);
    assert_eq!(g.apply_to(0.25), 0.5);
    assert_eq!(g.apply_to(0.5), 1.);
    assert_eq!(g.apply_to(0.75), 0.5);
    assert_eq!(g.apply_to(1.), 1.);
  }

  #[test]
  fn discrete_gradient() {
    let g = Gradient::Discrete {
      n: 2,
      gradient: Box::new(Gradient::Linear { factor: 1. }),
    };

    assert_eq!(g.apply_to(0.), 0.);
    assert_eq!(g.apply_to(0.25), 0.);
    assert_eq!(g.apply_to(0.5), 1.);
    assert_eq!(g.apply_to(0.75), 1.);
    assert_eq!(g.apply_to(1.), 1.);

    let g = Gradient::Discrete {
      n: 3,
      gradient: Box::new(Gradient::Linear { factor: 1. }),
    };

    assert_eq!(g.apply_to(0.), 0.);
    assert_eq!(g.apply_to(0.33), 0.);
    assert_eq!(g.apply_to(0.34), 0.5);
    assert_eq!(g.apply_to(0.66), 0.5);
    assert_eq!(g.apply_to(0.67), 1.);
    assert_eq!(g.apply_to(1.), 1.);

    let g = Gradient::Discrete {
      n: 4,
      gradient: Box::new(Gradient::Linear { factor: 1. }),
    };

    assert_eq!(g.apply_to(0.), 0.);
    assert_eq!(g.apply_to(0.24), 0.);
    assert_eq!((g.apply_to(0.25) * 100.).floor(), 33.);
    assert_eq!((g.apply_to(0.49) * 100.).floor(), 33.);
    assert_eq!((g.apply_to(0.5) * 100.).floor(), 66.);
    assert_eq!((g.apply_to(0.74) * 100.).floor(), 66.);
    assert_eq!(g.apply_to(0.75), 1.);
    assert_eq!(g.apply_to(1.), 1.);
  }

  #[test]
  fn grid_pos1() {
    let (x, y) = grid_pos(1.5, 0.5, 0., 2., 0., 2., 1., 1.).unwrap();
    assert_eq!(x, 1);
    assert_eq!(y, 0);
  }

  #[test]
  fn grid_pos2() {
    let p = grid_pos(1.5, -0.5, 0., 2., 0., 2., 1., 1.);
    assert_eq!(p, None);
  }

  #[test]
  fn grid_pos3() {
    let p = grid_pos(1.5, 2.5, 0., 2., 0., 2., 1., 1.);
    assert_eq!(p, None);
  }

  #[test]
  fn grid_pos4() {
    let (x, y) =
      grid_pos(0., 0., -50., 50., -50., 50., 1., 1.).unwrap();

    assert_eq!(x, 50);
    assert_eq!(y, 50);
  }

  #[test]
  fn grid_pos5() {
    let (x, y) =
      grid_pos(0., 0., -1., 1., -1., 1., 0.01, 0.005).unwrap();

    assert_eq!(x, 100);
    assert_eq!(y, 200);
  }

  #[test]
  fn grid_pos6() {
    let (x, y) =
      grid_pos(0.999, 0.999, -1., 1., -1., 1., 0.01, 0.005).unwrap();

    assert_eq!(x, 199);
    assert_eq!(y, 399);
  }
}
