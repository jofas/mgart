use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

use num::cast;

use crate::util::gradient::Gradient;

pub mod colors;

use colors::{Color, LCH, RGB};

#[derive(Serialize, Deserialize, DisplayAsJson, Clone)]
#[serde(from = "ColorMap1dDeserializer")]
pub struct ColorMap1d {
  map: Vec<LCH>,
  gradient: Gradient,
}

impl ColorMap1d {
  #[must_use]
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

  #[must_use]
  pub fn gradient(&self) -> &Gradient {
    &self.gradient
  }

  #[must_use]
  pub fn with_gradient(mut self, g: Gradient) -> Self {
    self.gradient = g;
    self
  }

  /// Creates an [`RGB`] color from `f`.
  ///
  /// # Panics
  ///
  /// Could panic, if something goes wrong with casting between
  /// [f64] and [usize].
  ///
  #[must_use]
  pub fn color(&self, f: f64) -> RGB {
    let f = self.gradient.apply(f);

    if f >= 1.0 {
      return self.map[self.map.len() - 1].rgb();
    }

    let interval = f * cast::<_, f64>(self.map.len() - 1).unwrap();

    let pos = interval.fract();

    let interval = cast::<_, usize>(interval).unwrap();

    let c1 = &self.map[interval];
    let c2 = &self.map[interval + 1];

    c1.interpolate(c2, pos).rgb()
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
