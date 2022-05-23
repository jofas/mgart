use cgmath::Vector3;

use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

use std::str::FromStr;

#[derive(Clone, Copy, Serialize, Deserialize, PartialEq, Debug)]
#[serde(from = "u32")]
#[serde(into = "u32")]
pub struct RgbaColor {
  r: u8,
  g: u8,
  b: u8,
  a: u8,
}

impl RgbaColor {
  /// Creates a new [RgbaColor] from hex.
  ///
  /// **Example:**
  ///
  /// ```rust
  /// use algorithmic_art::util::RgbaColor;
  ///
  /// // sunshine yellow
  /// let c = RgbaColor::new_hex(0xFFFD37FF);
  ///
  /// assert_eq!(c.r(), 255);
  /// assert_eq!(c.g(), 253);
  /// assert_eq!(c.b(), 55);
  /// assert_eq!(c.a(), 255);
  /// ```
  ///
  pub fn new_hex(color: u32) -> Self {
    Self {
      r: ((color & 0xFF000000) >> 24) as u8,
      g: ((color & 0xFF0000) >> 16) as u8,
      b: ((color & 0xFF00) >> 8) as u8,
      a: (color & 0xFF) as u8,
    }
  }

  /// Creates a new [RgbaColor] from rgba.
  ///
  /// **Example:**
  ///
  /// ```rust
  /// use algorithmic_art::util::RgbaColor;
  ///
  /// // sunshine yellow
  /// let c = RgbaColor::new_rgba(255, 253, 55, 255);
  ///
  /// assert_eq!(c.r(), 255);
  /// assert_eq!(c.g(), 253);
  /// assert_eq!(c.b(), 55);
  /// assert_eq!(c.a(), 255);
  ///
  /// assert_eq!(c.as_u32(), 0xFFFD37FF);
  /// ```
  ///
  pub fn new_rgba(r: u8, g: u8, b: u8, a: u8) -> Self {
    Self { r, g, b, a }
  }

  pub fn r(&self) -> u8 {
    self.r
  }

  pub fn g(&self) -> u8 {
    self.g
  }

  pub fn b(&self) -> u8 {
    self.b
  }

  pub fn a(&self) -> u8 {
    self.a
  }

  pub fn as_vec(&self) -> [u8; 4] {
    [self.r, self.g, self.b, self.a]
  }

  pub fn as_u32(&self) -> u32 {
    ((self.r() as u32) << 24)
      | ((self.g() as u32) << 16)
      | ((self.b() as u32) << 8)
      | (self.a() as u32)
  }
}

impl From<u32> for RgbaColor {
  fn from(x: u32) -> Self {
    Self::new_hex(x)
  }
}

impl Into<u32> for RgbaColor {
  fn into(self) -> u32 {
    self.as_u32()
  }
}

#[derive(Serialize, Deserialize, DisplayAsJson)]
pub struct ColorMap1D(Vec<RgbaColor>);

impl ColorMap1D {
  pub fn new(colors: Vec<RgbaColor>) -> Self {
    let colors = if colors.len() >= 2 {
      colors
    } else if colors.len() == 1 {
      vec![RgbaColor::new_hex(0xFFFFFFFF), colors[0]]
    } else {
      vec![
        RgbaColor::new_hex(0xFFFFFFFF),
        RgbaColor::new_hex(0x000000FF),
      ]
    };

    Self(colors)
  }

  pub fn value(&self, x: f32) -> RgbaColor {
    let x = 0.0_f32.max(0.99_f32.min(x));

    let interval = x * (self.0.len() - 1) as f32;
    let pos = interval.fract();

    let c1 = &self.0[interval as usize];
    let c2 = &self.0[interval as usize + 1];

    let v1: Vector3<f32> = self.color_to_vec3(&c1);
    let v2: Vector3<f32> = self.color_to_vec3(&c2);

    let res = (v2 - v1) * pos + v1;

    RgbaColor::new_rgba(
      res.x.abs() as u8,
      res.y.abs() as u8,
      res.z.abs() as u8,
      255,
    )
  }

  fn color_to_vec3(&self, c: &RgbaColor) -> Vector3<f32> {
    Vector3::new(c.r() as f32, c.g() as f32, c.b() as f32)
  }
}

impl FromStr for ColorMap1D {
  type Err = serde_json::Error;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    serde_json::from_str(s)
  }
}

#[cfg(test)]
mod tests {
  use super::{ColorMap1D, RgbaColor};

  #[test]
  fn deserialize_rgba_color() {
    let c: u32 = 0xFFFD37FF;

    let json = format!("{}", c);

    assert_eq!(
      serde_json::from_str::<RgbaColor>(&json).unwrap(),
      RgbaColor::new_hex(0xFFFD37FF),
    );
  }

  #[test]
  fn display_color_map_1d() {
    let cm = ColorMap1D::new(vec![]);

    assert_eq!(cm.to_string(), "[4294967295,255]");
  }
}
