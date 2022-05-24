use num_complex::Complex;

use clap::ArgEnum;

use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

use std::f64::consts::PI;
use std::fmt;
use std::str::FromStr;

pub mod colors;

use colors::{Color, LCH, RGB};

/// Representation of a complex number.
///
/// This is intended to be used as means for parsing user input,
/// not for doing calculations.
/// So [ComplexNumber] does not implement any math operations,
/// but supports the conversion to [Complex].
///
#[derive(Serialize, Deserialize)]
#[serde(untagged)]
pub enum ComplexNumber {
  Cartesian { re: f64, im: f64 },
  Polar { r: f64, theta: f64 },
}

impl FromStr for ComplexNumber {
  type Err = serde_json::Error;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    serde_json::from_str(s)
  }
}

impl Into<Complex<f64>> for &ComplexNumber {
  fn into(self) -> Complex<f64> {
    match self {
      ComplexNumber::Cartesian { re, im } => Complex::new(*re, *im),
      ComplexNumber::Polar { r, theta } => {
        Complex::from_polar(*r, *theta)
      }
    }
  }
}

/// How the [ColorMap1d] should be used.
///
/// Can be provided by the user as input and then passed to the
/// [ColorMap1d::color] method of a [color map](ColorMap1d).
///
#[derive(
  Serialize, Deserialize, ArgEnum, Clone, PartialEq, Debug,
)]
#[serde(rename_all = "lowercase")]
pub enum ColorMethod {
  #[clap(name = "linear")]
  Linear,
  #[clap(name = "sine")]
  Sine,
}

impl FromStr for ColorMethod {
  type Err = serde_json::Error;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    serde_json::from_str(&format!("\"{}\"", s))
  }
}

impl fmt::Display for ColorMethod {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self {
      Self::Linear => write!(f, "linear"),
      Self::Sine => write!(f, "sine"),
    }
  }
}

#[derive(Serialize, Deserialize, DisplayAsJson)]
#[serde(from = "ColorMap1dDeserializer")]
pub struct ColorMap1d {
  map: Vec<LCH>,
  method: ColorMethod,
}

impl ColorMap1d {
  pub fn new(map: Vec<Color>, method: ColorMethod) -> Self {
    let map = if map.len() >= 2 {
      map
    } else if map.len() == 1 {
      vec![Color::WHITE, map[0]]
    } else {
      vec![Color::WHITE, Color::BLACK]
    };

    let map: Vec<LCH> = map.into_iter().map(|c| c.lch()).collect();

    Self { map, method }
  }

  pub fn color(&self, x: f64) -> RGB {
    match self.method {
      ColorMethod::Linear => self.linear(x),
      ColorMethod::Sine => self.sine(x),
    }
  }

  fn linear(&self, f: f64) -> RGB {
    let f = f.clamp(0., 1.);

    if 1.0 - f <= f64::EPSILON {
      return self.map[self.map.len() - 1].rgb();
    }

    let interval = f * (self.map.len() - 1) as f64;
    let pos = interval.fract() as f64;

    let c1 = &self.map[interval as usize];
    let c2 = &self.map[interval as usize + 1];

    c1.interpolate(&c2, pos).rgb()
  }

  /// Computes the color at point `x` by passing `sin(x * PI)` into
  /// [Self::linear].
  ///
  fn sine(&self, f: f64) -> RGB {
    let f = f.clamp(0., 1.);
    self.linear((f * 2. * PI).sin().abs())
  }
}

impl FromStr for ColorMap1d {
  type Err = serde_json::Error;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    serde_json::from_str(s)
  }
}

#[derive(Deserialize)]
struct ColorMap1dDeserializer {
  map: Vec<Color>,
  method: ColorMethod,
}

impl From<ColorMap1dDeserializer> for ColorMap1d {
  fn from(cm: ColorMap1dDeserializer) -> Self {
    Self::new(cm.map, cm.method)
  }
}
