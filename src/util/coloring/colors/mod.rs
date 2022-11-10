use serde::{Deserialize, Serialize};

use std::f64::consts::PI;

mod serde_nan;

#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
#[serde(rename_all = "lowercase")]
#[serde(tag = "type")]
pub enum Color {
  LCH(LCH),
  RGB(RGB),
}

impl Color {
  pub const BLACK: Self = Self::LCH(LCH::BLACK);
  pub const WHITE: Self = Self::LCH(LCH::WHITE);

  pub fn lch(&self) -> LCH {
    match self {
      Self::LCH(lch) => *lch,
      Self::RGB(rgb) => rgb.lch(),
    }
  }

  pub fn rgb(&self) -> RGB {
    match self {
      Self::LCH(lch) => lch.rgb(),
      Self::RGB(rgb) => *rgb,
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

  /// Returns the [RGB] representation of the color defined by
  /// `self`.
  ///
  pub fn rgb(&self) -> RGB {
    self.lab().rgb()
  }

  /// Returns the [LAB] representation of the color defined by `self`.
  ///
  fn lab(&self) -> LAB {
    let h = if self.h().is_nan() { 0_f64 } else { self.h() };
    let h = h * PI / 180.;

    LAB::new(self.l(), h.cos() * self.c(), h.sin() * self.c())
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

  /// Rounds the [l](Self::l), [c](Self::c) and [h](Self::h) values
  /// to the nearest number with the provided precision.
  ///
  #[allow(dead_code)]
  fn round(&self, digits: i32) -> Self {
    let pow = 10_f64.powi(digits);

    let l = (self.l() * pow).round() / pow;
    let c = (self.c() * pow).round() / pow;
    let h = (self.h() * pow).round() / pow;

    Self::new(l, c, h)
  }
}

impl PartialEq for LCH {
  fn eq(&self, other: &Self) -> bool {
    let lc_eq = self.l() == other.l() && self.c() == other.c();

    if self.h() != other.h() {
      self.h().is_nan() && other.h().is_nan() && lc_eq
    } else {
      lc_eq
    }
  }
}

#[derive(
  Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Debug,
)]
pub struct RGB {
  r: u8,
  g: u8,
  b: u8,
}

impl RGB {
  /// Creates a new [RGB] color.
  ///
  pub const fn new(r: u8, g: u8, b: u8) -> Self {
    Self { r, g, b }
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

  pub fn lch(&self) -> LCH {
    self.lab().lch()
  }

  pub fn as_vec(&self) -> [u8; 3] {
    [self.r, self.g, self.b]
  }

  fn lab(&self) -> LAB {
    let r = Self::rgb_xyz(self.r());
    let g = Self::rgb_xyz(self.g());
    let b = Self::rgb_xyz(self.b());

    let x = Self::xyz_lab(
      (0.4124564 * r + 0.3575761 * g + 0.1804375 * b) / LAB::XN,
    );
    let y = Self::xyz_lab(
      (0.2126729 * r + 0.7151522 * g + 0.0721750 * b) / LAB::YN,
    );
    let z = Self::xyz_lab(
      (0.0193339 * r + 0.1191920 * g + 0.9503041 * b) / LAB::ZN,
    );

    let l = (116. * y - 16.).max(0.);

    LAB::new(l, 500. * (x - y), 200. * (y - z))
  }

  fn rgb_xyz(r: u8) -> f64 {
    let r = r as f64 / 255.;

    if r <= 0.04045 {
      r / 12.92
    } else {
      ((r + 0.055) / 1.055).powf(2.4)
    }
  }

  fn xyz_lab(t: f64) -> f64 {
    if t > LAB::T3 {
      t.powf(1. / 3.)
    } else {
      t / LAB::T2 + LAB::T0
    }
  }
}

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, PartialEq)]
struct LAB {
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
  const T3: f64 = 0.008856452;

  const fn new(l: f64, a: f64, b: f64) -> Self {
    Self { l, a, b }
  }

  /// The lightness of the color.
  ///
  fn l(&self) -> f64 {
    self.l
  }

  /// The value on the a-axis of the color.
  ///
  fn a(&self) -> f64 {
    self.a
  }

  /// The value on the b-axis of the color.
  ///
  fn b(&self) -> f64 {
    self.b
  }

  /// Rounds the [l](Self::l), [a](Self::a) and [b](Self::b) values
  /// to the nearest number with the provided precision.
  ///
  #[allow(dead_code)]
  fn round(&self, digits: i32) -> Self {
    let pow = 10_f64.powi(digits);

    let l = (self.l() * pow).round() / pow;
    let a = (self.a() * pow).round() / pow;
    let b = (self.b() * pow).round() / pow;

    Self::new(l, a, b)
  }

  /// Returns the [LCH] representation of the color defined by
  /// `self`.
  ///
  fn lch(&self) -> LCH {
    let c = (self.a().powi(2) + self.b().powi(2)).sqrt();

    let h = if 0. + (c * 1.0e4).round() <= f64::EPSILON {
      f64::NAN
    } else {
      (self.b().atan2(self.a()) * 180. / PI + 360.) % 360.
    };

    LCH::new(self.l(), c, h)
  }

  /// Returns the [RGB] representation of the color defined by
  /// `self`.
  ///
  fn rgb(&self) -> RGB {
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

    RGB::new(r, g, b)
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

#[cfg(test)]
mod test {
  use super::{LAB, LCH, RGB};

  use std::f64;

  const BLACK_RGB: RGB = RGB::new(0, 0, 0);
  const WHITE_RGB: RGB = RGB::new(255, 255, 255);
  const GREY_RGB: RGB = RGB::new(128, 128, 128);
  const RED_RGB: RGB = RGB::new(255, 0, 0);
  const YELLOW_RGB: RGB = RGB::new(255, 255, 0);
  const GREEN_RGB: RGB = RGB::new(0, 255, 0);
  const CYAN_RGB: RGB = RGB::new(0, 255, 255);
  const BLUE_RGB: RGB = RGB::new(0, 0, 255);
  const MAGENTA_RGB: RGB = RGB::new(255, 0, 255);

  const BLACK_LCH: LCH = LCH::new(0., 0., f64::NAN);
  const WHITE_LCH: LCH = LCH::new(100., 0., f64::NAN);
  const GREY_LCH: LCH = LCH::new(53.59, 0., f64::NAN);
  const RED_LCH: LCH = LCH::new(53.24, 104.55, 40.);
  const YELLOW_LCH: LCH = LCH::new(97.14, 96.91, 102.85);
  const GREEN_LCH: LCH = LCH::new(87.73, 119.78, 136.02);
  const CYAN_LCH: LCH = LCH::new(91.11, 50.12, 196.37);
  const BLUE_LCH: LCH = LCH::new(32.3, 133.81, 306.28);
  const MAGENTA_LCH: LCH = LCH::new(60.32, 115.54, 328.23);

  const BLACK_LAB: LAB = LAB::new(0., 0., 0.);
  const WHITE_LAB: LAB = LAB::new(100., 0., 0.);
  const GREY_LAB: LAB = LAB::new(53.59, 0., 0.);
  const RED_LAB: LAB = LAB::new(53.24, 80.09, 67.2);
  const YELLOW_LAB: LAB = LAB::new(97.14, -21.55, 94.48);
  const GREEN_LAB: LAB = LAB::new(87.73, -86.19, 83.18);
  const CYAN_LAB: LAB = LAB::new(91.11, -48.09, -14.13);
  const BLUE_LAB: LAB = LAB::new(32.3, 79.18, -107.87);
  const MAGENTA_LAB: LAB = LAB::new(60.32, 98.23, -60.83);

  #[test]
  fn lch2lab() {
    assert_eq!(BLACK_LCH.lab().round(2), BLACK_LAB);
    assert_eq!(WHITE_LCH.lab().round(2), WHITE_LAB);
    assert_eq!(GREY_LCH.lab().round(2), GREY_LAB);
    assert_eq!(RED_LCH.lab().round(2), RED_LAB);
    assert_eq!(YELLOW_LCH.lab().round(2), YELLOW_LAB);
    assert_eq!(GREEN_LCH.lab().round(2), GREEN_LAB);
    assert_eq!(CYAN_LCH.lab().round(2), CYAN_LAB);
    assert_eq!(BLUE_LCH.lab().round(2), BLUE_LAB);
    assert_eq!(MAGENTA_LCH.lab().round(2), MAGENTA_LAB);
  }

  #[test]
  fn lab2lch() {
    assert_eq!(BLACK_LAB.lch().round(2), BLACK_LCH);
    assert_eq!(WHITE_LAB.lch().round(2), WHITE_LCH);
    assert_eq!(GREY_LAB.lch().round(2), GREY_LCH);
    assert_eq!(RED_LAB.lch().round(2), RED_LCH);
    assert_eq!(YELLOW_LAB.lch().round(2), YELLOW_LCH);
    assert_eq!(GREEN_LAB.lch().round(2), GREEN_LCH);
    assert_eq!(CYAN_LAB.lch().round(2), CYAN_LCH);
    assert_eq!(BLUE_LAB.lch().round(2), BLUE_LCH);
    assert_eq!(MAGENTA_LAB.lch().round(2), MAGENTA_LCH);
  }

  #[test]
  fn lab2rgb() {
    assert_eq!(BLACK_LAB.rgb(), BLACK_RGB);
    assert_eq!(WHITE_LAB.rgb(), WHITE_RGB);
    assert_eq!(GREY_LAB.rgb(), GREY_RGB);
    assert_eq!(RED_LAB.rgb(), RED_RGB);
    assert_eq!(YELLOW_LAB.rgb(), YELLOW_RGB);
    assert_eq!(GREEN_LAB.rgb(), GREEN_RGB);
    assert_eq!(CYAN_LAB.rgb(), CYAN_RGB);
    assert_eq!(BLUE_LAB.rgb(), BLUE_RGB);
    assert_eq!(MAGENTA_LAB.rgb(), MAGENTA_RGB);
  }

  #[test]
  fn lch2rgb() {
    assert_eq!(BLACK_LCH.rgb(), BLACK_RGB);
    assert_eq!(WHITE_LCH.rgb(), WHITE_RGB);
    assert_eq!(GREY_LCH.rgb(), GREY_RGB);
    assert_eq!(RED_LCH.rgb(), RED_RGB);
    assert_eq!(YELLOW_LCH.rgb(), YELLOW_RGB);
    assert_eq!(GREEN_LCH.rgb(), GREEN_RGB);
    assert_eq!(CYAN_LCH.rgb(), CYAN_RGB);
    assert_eq!(BLUE_LCH.rgb(), BLUE_RGB);
    assert_eq!(MAGENTA_LCH.rgb(), MAGENTA_RGB);
  }

  #[test]
  fn rgb2lab() {
    assert_eq!(BLACK_RGB.lab().round(2), BLACK_LAB);
    assert_eq!(WHITE_RGB.lab().round(2), WHITE_LAB);
    assert_eq!(GREY_RGB.lab().round(2), GREY_LAB);
    assert_eq!(RED_RGB.lab().round(2), RED_LAB);
    assert_eq!(YELLOW_RGB.lab().round(2), YELLOW_LAB);
    assert_eq!(
      GREEN_RGB.lab().round(2),
      LAB::new(87.73, -86.18, 83.18)
    );
    assert_eq!(CYAN_RGB.lab().round(2), CYAN_LAB);
    assert_eq!(
      BLUE_RGB.lab().round(2),
      LAB::new(32.3, 79.19, -107.86)
    );
    assert_eq!(
      MAGENTA_RGB.lab().round(2),
      LAB::new(60.32, 98.23, -60.82)
    );
  }

  #[test]
  fn rgb2lch() {
    assert_eq!(BLACK_RGB.lch().round(2), BLACK_LCH);
    assert_eq!(WHITE_RGB.lch().round(2), WHITE_LCH);
    assert_eq!(GREY_RGB.lch().round(2), GREY_LCH);
    assert_eq!(RED_RGB.lch().round(2), RED_LCH);
    assert_eq!(YELLOW_RGB.lch().round(2), YELLOW_LCH);
    assert_eq!(GREEN_RGB.lch().round(2), GREEN_LCH);
    assert_eq!(
      CYAN_RGB.lch().round(2),
      LCH::new(91.11, 50.12, 196.38)
    );
    assert_eq!(BLUE_RGB.lch().round(2), BLUE_LCH);
    assert_eq!(MAGENTA_RGB.lch().round(2), MAGENTA_LCH);
  }
}
