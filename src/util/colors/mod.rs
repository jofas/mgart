use serde::{Deserialize, Serialize};

use std::f64::consts::PI;

mod serde_nan;

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "lowercase")]
pub enum ColorSpace {
  LCH,
  RGBA,
}

#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
#[serde(rename_all = "lowercase")]
#[serde(tag = "type")]
pub enum Color {
  LCH(LCH),
  RGBA(RGBA),
}

impl Color {
  pub const BLACK: Self = Self::LCH(LCH::BLACK);
  pub const WHITE: Self = Self::LCH(LCH::WHITE);

  pub fn as_rgba(&self) -> Self {
    match self {
      Self::LCH(lch) => Self::RGBA(lch.rgba()),
      Self::RGBA(_) => *self,
    }
  }

  pub fn rgba(&self) -> RGBA {
    match self {
      Self::LCH(lch) => lch.rgba(),
      Self::RGBA(rgba) => *rgba,
    }
  }

  pub fn as_lch(&self) -> Self {
    match self {
      Self::LCH(_) => *self,
      Self::RGBA(rgba) => Self::LCH(rgba.lch()),
    }
  }

  pub fn lch(&self) -> LCH {
    match self {
      Self::LCH(lch) => *lch,
      Self::RGBA(rgba) => rgba.lch(),
    }
  }

  /// Interpolates `self` with `other`.
  ///
  /// The interpolation happens in the color space of `self`.
  ///
  pub fn interpolate(&self, other: &Self, f: f64) -> Color {
    match self {
      Self::LCH(lch) => Self::LCH(lch.interpolate(&other.lch(), f)),
      Self::RGBA(rgba) => {
        Self::RGBA(rgba.interpolate(&other.rgba(), f))
      }
    }
  }
}

#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
pub struct LCH {
  #[serde(with = "serde_nan")]
  l: f64,
  #[serde(with = "serde_nan")]
  c: f64,
  #[serde(with = "serde_nan")]
  h: f64,
}

impl LCH {
  pub const BLACK: LCH = Self::new(0., 0., f64::NAN);
  pub const WHITE: LCH = Self::new(100., 0., f64::NAN);

  pub const fn new(l: f64, c: f64, h: f64) -> Self {
    Self { l, c, h }
  }

  /// The lightness of the color.
  ///
  pub fn l(&self) -> f64 {
    self.l
  }

  /// The chroma of the color.
  ///
  pub fn c(&self) -> f64 {
    self.c
  }

  /// The hue of the color.
  ///
  pub fn h(&self) -> f64 {
    self.h
  }

  /// Returns the [LAB] representation of the color defined by `self`.
  ///
  pub fn lab(&self) -> LAB {
    let h = if self.h().is_nan() { 0_f64 } else { self.h() };
    let h = h * PI / 180.;

    return LAB::new(
      self.l(),
      h.cos() * self.c(),
      h.sin() * self.c(),
    );
  }

  /// Returns the [RGBA] representation of the color defined by
  /// `self`.
  ///
  pub fn rgba(&self) -> RGBA {
    self.lab().rgba()
  }

  pub fn interpolate(&self, other: &LCH, f: f64) -> LCH {
    let f = f.clamp(0., 1.);

    let h0 = self.h();
    let h1 = other.h();

    let c0 = self.c();
    let c1 = other.c();

    let l0 = self.l();
    let l1 = other.l();

    let mut h = f64::NAN;
    let mut c = None;

    if !h0.is_nan() && !h1.is_nan() {
      let dh = if h1 > h0 && h1 - h0 > 180. {
        h1 - h0 + 360.
      } else if h1 < h0 && h0 - h1 > 180. {
        h1 + 360. - h0
      } else {
        h1 - h0
      };

      h = h0 + f * dh;
    } else if !h0.is_nan() {
      h = h0;
      if l1 == 1. || l1 == 0. {
        c = Some(c0);
      }
    } else if !h1.is_nan() {
      h = h1;
      if l0 == 1. || l0 == 0. {
        c = Some(c1);
      }
    }

    let c = match c {
      Some(c) => c,
      None => c0 + f * (c1 - c0),
    };

    let l = l0 + f * (l1 - l0);

    LCH::new(l, c, h)
  }
}

#[derive(Debug, PartialEq)]
pub struct LAB {
  l: f64,
  a: f64,
  b: f64,
}

impl LAB {
  const XN: f64 = 0.950470;
  const YN: f64 = 1.;
  const ZN: f64 = 1.088830;

  const T0: f64 = 0.137931034;
  const T1: f64 = 0.206896552;
  const T2: f64 = 0.12841855;

  pub fn new(l: f64, a: f64, b: f64) -> Self {
    Self { l, a, b }
  }

  /// The lightness of the color.
  ///
  pub fn l(&self) -> f64 {
    self.l
  }

  /// The value on the a-axis of the color.
  ///
  pub fn a(&self) -> f64 {
    self.a
  }

  /// The value on the b-axis of the color.
  ///
  pub fn b(&self) -> f64 {
    self.b
  }

  /// Rounds the [l](Self::l), [a](Self::a) and [b](Self::b) values
  /// to the nearest number with the provided precision.
  ///
  /// **Example:**
  ///
  /// ```rust
  /// use algorithmic_art::util::colors::LAB;
  ///
  /// let c = LAB::new(0.005, 0.0049, 1.2);
  ///
  /// assert_eq!(c.round(2), LAB::new(0.01, 0., 1.2));
  /// ```
  ///
  pub fn round(&self, digits: i32) -> Self {
    let pow = 10_f64.powi(digits);

    let l = (self.l() * pow).round() / pow;
    let a = (self.a() * pow).round() / pow;
    let b = (self.b() * pow).round() / pow;

    Self::new(l, a, b)
  }

  /// Returns the [RGBA] representation of the color defined by
  /// `self`.
  ///
  pub fn rgba(&self) -> RGBA {
    let y = (self.l() + 16.) / 116.;

    let x = if self.a().is_nan() {
      y
    } else {
      y + self.a() / 500.
    };

    let z = if self.a().is_nan() {
      y
    } else {
      y - self.b() / 200.
    };

    let y = Self::YN * Self::lab_xyz(y);
    let x = Self::XN * Self::lab_xyz(x);
    let z = Self::ZN * Self::lab_xyz(z);

    let r =
      Self::xyz_rgb(3.2404542 * x - 1.5371385 * y - 0.4985314 * z);
    let g =
      Self::xyz_rgb(-0.9692660 * x + 1.8760108 * y + 0.0415560 * z);
    let b =
      Self::xyz_rgb(0.0556434 * x - 0.2040259 * y + 1.0572252 * z);

    RGBA::new_rgba(r, g, b, 255)
  }

  fn xyz_rgb(r: f64) -> u8 {
    let x = if r <= 0.00304 {
      12.92 * r
    } else {
      1.055 * r.powf(1. / 2.4) - 0.055
    };

    (255. * x.clamp(0., 1.)).round() as u8
  }

  fn lab_xyz(t: f64) -> f64 {
    if t > Self::T1 {
      t.powi(3)
    } else {
      Self::T2 * (t - Self::T0)
    }
  }
}

#[derive(Clone, Copy, Serialize, Deserialize, PartialEq, Debug)]
pub struct RGBA {
  r: u8,
  g: u8,
  b: u8,
  a: u8,
}

impl RGBA {
  /// Creates a new [RGBA] from hex.
  ///
  /// **Example:**
  ///
  /// ```rust
  /// use algorithmic_art::util::colors::RGBA;
  ///
  /// // sunshine yellow
  /// let c = RGBA::new_hex(0xFFFD37FF);
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

  /// Creates a new [RGBA] from rgba.
  ///
  /// **Example:**
  ///
  /// ```rust
  /// use algorithmic_art::util::colors::RGBA;
  ///
  /// // sunshine yellow
  /// let c = RGBA::new_rgba(255, 253, 55, 255);
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

  pub fn lch(&self) -> LCH {
    // TODO
    unimplemented!()
  }

  pub fn interpolate(&self, other: &Self, f: f64) -> RGBA {
    let f = f.clamp(0., 1.);

    let r0 = self.r() as f64;
    let r1 = other.r() as f64;

    let g0 = self.g() as f64;
    let g1 = other.g() as f64;

    let b0 = self.b() as f64;
    let b1 = other.b() as f64;

    let a0 = self.a() as f64;
    let a1 = other.a() as f64;

    let r = r0 + f * (r1 - r0);
    let g = g0 + f * (g1 - g0);
    let b = b0 + f * (b1 - b0);
    let a = a0 + f * (a1 - a0);

    Self::new_rgba(
      r.abs() as u8,
      g.abs() as u8,
      b.abs() as u8,
      a.abs() as u8,
    )
  }

  pub fn as_vec(&self) -> [u8; 4] {
    [self.r, self.g, self.b, self.a]
  }
}

#[cfg(test)]
mod test {
  use super::{LAB, LCH, RGBA};

  use std::f64;

  #[test]
  fn lch2lab() {
    let black = LCH::new(0., 0., f64::NAN);
    let white = LCH::new(100., 0., f64::NAN);
    let grey = LCH::new(53.59, 0., f64::NAN);
    let red = LCH::new(53.24, 104.55, 40.);
    let yellow = LCH::new(97.14, 96.91, 102.85);
    let green = LCH::new(87.73, 119.78, 136.02);
    let cyan = LCH::new(91.11, 50.12, 196.38);
    let blue = LCH::new(32.3, 133.81, 306.28);
    let magenta = LCH::new(60.32, 115.54, 328.23);

    assert_eq!(black.lab().round(2), LAB::new(0., 0., 0.));
    assert_eq!(white.lab().round(2), LAB::new(100., 0., 0.));
    assert_eq!(grey.lab().round(2), LAB::new(53.59, 0., 0.));
    assert_eq!(red.lab().round(2), LAB::new(53.24, 80.09, 67.2));
    assert_eq!(yellow.lab().round(2), LAB::new(97.14, -21.55, 94.48));
    assert_eq!(green.lab().round(2), LAB::new(87.73, -86.19, 83.18));
    assert_eq!(cyan.lab().round(2), LAB::new(91.11, -48.09, -14.13));
    assert_eq!(blue.lab().round(2), LAB::new(32.3, 79.18, -107.87));
    assert_eq!(
      magenta.lab().round(2),
      LAB::new(60.32, 98.23, -60.83)
    );
  }

  #[test]
  fn lab2rgba() {
    let black = LAB::new(0., 0., 0.);
    let white = LAB::new(100., 0., 0.);
    let grey = LAB::new(53.59, 0., 0.);
    let red = LAB::new(53.24, 80.09, 67.2);
    let yellow = LAB::new(97.14, -21.55, 94.48);
    let green = LAB::new(87.73, -86.19, 83.18);
    let cyan = LAB::new(91.11, -48.09, -14.13);
    let blue = LAB::new(32.3, 79.18, -107.87);
    let magenta = LAB::new(60.32, 98.23, -60.83);

    assert_eq!(black.rgba(), RGBA::new_rgba(0, 0, 0, 255));
    assert_eq!(white.rgba(), RGBA::new_rgba(255, 255, 255, 255));
    assert_eq!(grey.rgba(), RGBA::new_rgba(128, 128, 128, 255));
    assert_eq!(red.rgba(), RGBA::new_rgba(255, 0, 0, 255));
    assert_eq!(yellow.rgba(), RGBA::new_rgba(255, 255, 0, 255));
    assert_eq!(green.rgba(), RGBA::new_rgba(0, 255, 0, 255));
    assert_eq!(cyan.rgba(), RGBA::new_rgba(0, 255, 255, 255));
    assert_eq!(blue.rgba(), RGBA::new_rgba(0, 0, 255, 255));
    assert_eq!(magenta.rgba(), RGBA::new_rgba(255, 0, 255, 255));
  }

  #[test]
  fn lch2rgba() {
    let black = LCH::new(0., 0., f64::NAN);
    let white = LCH::new(100., 0., f64::NAN);
    let grey = LCH::new(53.59, 0., f64::NAN);
    let red = LCH::new(53.24, 104.55, 40.);
    let yellow = LCH::new(97.14, 96.91, 102.85);
    let green = LCH::new(87.73, 119.78, 136.02);
    let cyan = LCH::new(91.11, 50.12, 196.38);
    let blue = LCH::new(32.3, 133.81, 306.28);
    let magenta = LCH::new(60.32, 115.54, 328.23);

    assert_eq!(black.rgba(), RGBA::new_rgba(0, 0, 0, 255));
    assert_eq!(white.rgba(), RGBA::new_rgba(255, 255, 255, 255));
    assert_eq!(grey.rgba(), RGBA::new_rgba(128, 128, 128, 255));
    assert_eq!(red.rgba(), RGBA::new_rgba(255, 0, 0, 255));
    assert_eq!(yellow.rgba(), RGBA::new_rgba(255, 255, 0, 255));
    assert_eq!(green.rgba(), RGBA::new_rgba(0, 255, 0, 255));
    assert_eq!(cyan.rgba(), RGBA::new_rgba(0, 255, 255, 255));
    assert_eq!(blue.rgba(), RGBA::new_rgba(0, 0, 255, 255));
    assert_eq!(magenta.rgba(), RGBA::new_rgba(255, 0, 255, 255));
  }
}
