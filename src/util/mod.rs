use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

use num_complex::Complex64;

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
      ComplexNumber::Polar { r, theta } => Complex64::from_polar(*r, *theta),
    }
  }
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "type")]
pub enum Gradient {
  Linear { factor: f64 }, // TODO: option to disable repeat
  Sin { factor: f64 },
  Inverted { gradient: Box<Gradient> },
  Wave { factor: f64 },
  Exp { exponent: f64 },
  SinExp { factor: f64 },
  Log { factor: f64 },
  Tanh { factor: f64 },
  // TODO: smoothstep with order, sine-ramp, discrete, b-spline
}

impl Gradient {
  pub fn apply_to(&self, f: f64) -> f64 {
    match self {
      Self::Linear { factor } => (f * factor).fract(),
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
      Self::Log { factor } => (f * factor + 1.).ln() / (factor + 1.).ln(),
      Self::Tanh { factor } => (f * factor).tanh() / factor.tanh(),
    }
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

    if 1.0 - f <= f64::EPSILON {
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
