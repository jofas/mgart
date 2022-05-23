use cgmath::Vector3;

use serde::Deserialize;

use std::str::FromStr;

#[derive(Clone, Copy, Deserialize)]
#[serde(from = "u32")]
pub struct RgbaColor {
  r: u8,
  g: u8,
  b: u8,
  a: u8,
}

impl RgbaColor {
  pub fn new_hex(color: u32) -> Self {
    Self {
      r: ((color & 0xFF000000) >> 24) as u8,
      g: ((color & 0xFF0000) >> 16) as u8,
      b: ((color & 0xFF00) >> 8) as u8,
      a: (color & 0xFF) as u8,
    }
  }

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
}

impl From<u32> for RgbaColor {
  fn from(x: u32) -> Self {
    Self::new_hex(x)
  }
}

#[derive(Deserialize)]
pub struct ColorMap1D {
  #[serde(flatten)]
  colors: Vec<RgbaColor>,
}

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

    Self { colors }
  }

  pub fn value(&self, x: f32) -> RgbaColor {
    let x = 0.0_f32.max(0.99_f32.min(x));

    let interval = x * (self.colors.len() - 1) as f32;
    let pos = interval.fract();

    let c1 = &self.colors[interval as usize];
    let c2 = &self.colors[interval as usize + 1];

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
