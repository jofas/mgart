use rayon::iter::{
  IndexedParallelIterator, IntoParallelIterator,
  IntoParallelRefMutIterator, ParallelIterator,
};
use rayon::slice::ParallelSliceMut;

use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

use num_complex::Complex64;

use map_macro::vec_no_clone;

use rand::random;
use rand_distr::{Distribution, Normal};

use std::f64::consts::PI;
use std::sync::atomic::{AtomicU32, Ordering};

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

impl Default for Gradient {
  fn default() -> Self {
    Self::Linear { factor: 1. }
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

  pub fn gradient(&self) -> &Gradient {
    &self.gradient
  }

  pub fn with_gradient(mut self, g: Gradient) -> Self {
    self.gradient = g;
    self
  }

  pub fn color(&self, f: f64) -> RGB {
    let f = self.gradient.apply_to(f);

    if f >= 1.0 {
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

#[derive(Serialize, Deserialize, DisplayAsJson)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "type")]
pub enum Smoothing {
  NonLocalMeans { n: usize, window_size: usize },
}

impl Smoothing {
  pub fn smooth(
    &self,
    buffer: &[f64],
    variance: &[f64],
    width: usize,
    height: usize,
  ) -> Vec<f64> {
    match self {
      Self::NonLocalMeans { n, window_size } => {
        let num_pixel = buffer.len();

        let sat = summed_area_table(&buffer, width);

        let mut smoothed: Vec<f64> = vec![0.; num_pixel];

        let processed_pixels = AtomicU32::new(0);

        smoothed.par_iter_mut().enumerate().for_each(|(i, pixel)| {
          let x = i % width;
          let y = i / width;

          let (wx0, wy0, wx1, wy1) = discrete_bounded_square(
            x,
            y,
            *window_size,
            width - 1,
            height - 1,
          );

          let (x0, y0, x1, y1) =
            discrete_bounded_square(x, y, *n, width - 1, height - 1);

          let c00 = sat[y0 * width + x0];
          let c01 = sat[y0 * width + x1];
          let c10 = sat[y1 * width + x0];
          let c11 = sat[y1 * width + x1];

          let bp = c00 + c11 - c01 - c10;

          if bp.is_nan() {
            dbg!(x, y, x0, y0, x1, y1, c00, c01, c10, c11);
            panic!("bp is nan");
          }

          let bp = bp / (x1 - x0 + 1) as f64 / (y1 - y0 + 1) as f64;

          let mut s = 0.;
          let mut cp = 0.;

          for x in wx0..=wx1 {
            for y in wy0..=wy1 {
              let (x0, y0, x1, y1) = discrete_bounded_square(
                x,
                y,
                *n,
                width - 1,
                height - 1,
              );

              let bq = sat[y1 * width + x1] + sat[y0 * width + x0]
                - sat[y1 * width + x0]
                - sat[y0 * width + x1];
              let bq =
                bq / (x1 - x0 + 1) as f64 / (y1 - y0 + 1) as f64;

              let fpq = (-((bq - bp).powi(2) / variance[i])).exp();

              s += buffer[i] * fpq;
              cp += fpq;
            }
          }

          *pixel = s / cp;

          let pc = processed_pixels.fetch_add(1, Ordering::SeqCst);
          print_progress(pc, num_pixel as u32, 100);
        });

        smoothed
      }
    }
  }
}

pub trait Sampling {
  type Output;

  fn sample(&self) -> Self::Output;
}

pub struct KDE<T, K> {
  elems: Vec<T>,
  probabilities: Vec<f64>,
  kernel: K,
}

impl<T, K: Fn(&T) -> T> KDE<T, K> {
  pub fn new(samples: Vec<(T, f64)>, kernel: K) -> Self {
    let (elems, probabilities) = samples.into_iter().unzip();

    Self {
      elems: elems,
      probabilities: probabilities,
      kernel,
    }
  }

  pub fn average_probability(&self) -> f64 {
    self.probabilities.iter().sum::<f64>() / self.elems.len() as f64
  }
}

impl<T, K: Fn(&T) -> T> Sampling for KDE<T, K> {
  type Output = T;

  fn sample(&self) -> Self::Output {
    loop {
      let idx = (random::<f64>() * self.elems.len() as f64) as usize;

      if random::<f64>() <= self.probabilities[idx] {
        return (self.kernel)(&self.elems[idx]);
      }
    }
  }
}

pub struct Viewport {
  pub x_min: f64,
  pub y_min: f64,
  pub x_max: f64,
  pub y_max: f64,
}

impl Viewport {
  pub fn new(x_min: f64, y_min: f64, x_max: f64, y_max: f64) -> Self {
    Self {
      x_min,
      y_min,
      x_max,
      y_max,
    }
  }

  pub fn contains_point(&self, x: f64, y: f64) -> bool {
    self.x_min <= x
      && x < self.x_max
      && self.y_min <= y
      && y < self.y_max
  }
}

pub fn grid_pos(
  x: f64,
  y: f64,
  delta_x: f64,
  delta_y: f64,
  viewport: &Viewport,
) -> Option<(usize, usize)> {
  if viewport.contains_point(x, y) {
    let x = ((x - viewport.x_min) / delta_x) as usize;
    let y = ((y - viewport.y_min) / delta_y) as usize;

    Some((x, y))
  } else {
    None
  }
}

pub fn discrete_bounded_square(
  x: usize,
  y: usize,
  n: usize,
  ubx: usize,
  uby: usize,
) -> (usize, usize, usize, usize) {
  let nh = n / 2;

  let xupper = (x + nh).min(ubx);
  let yupper = (y + nh).min(uby);

  let (x, y) = if n % 2 == 0 { (x + 1, y + 1) } else { (x, y) };

  let xlower = if nh > x + 1 { 0 } else { x + 1 - nh };
  let ylower = if nh > y + 1 { 0 } else { y + 1 - nh };

  (xlower, ylower, xupper, yupper)
}

pub fn random_complex() -> Complex64 {
  Complex64::from_polar(
    2. * random::<f64>(),
    2. * PI * random::<f64>(),
  )
}

pub fn min_max(v: &[f64]) -> (f64, f64) {
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

pub fn summed_area_table(grid: &[f64], width: usize) -> Vec<f64> {
  let mut sat = grid.to_vec();

  for i in 1..sat.len() {
    let x = i % width;
    let y = i / width;

    if x > 0 {
      sat[i] += sat[i - 1];
    }

    if y > 0 {
      sat[i] += sat[(y - 1) * width + x];
    }

    if x > 0 && y > 0 {
      sat[i] -= sat[(y - 1) * width + x - 1];
    }
  }

  sat
}

pub fn print_progress(i: u32, n: u32, interval: u32) {
  if i % interval == interval - 1 || i == n - 1 {
    let p = i as f64 / n as f64 * 100.;
    print!("{}/{} iterations done ({:.2}%)\r", i + 1, n, p);
  }
}

#[cfg(test)]
mod tests {
  use super::{grid_pos, summed_area_table, Gradient};

  #[test]
  fn sat() {
    let image = [1.; 16];
    let width = 4;

    let sat = summed_area_table(&image, width);

    assert_eq!(
      sat,
      vec![
        1., 2., 3., 4., 2., 4., 6., 8., 3., 6., 9., 12., 4., 8., 12.,
        16.,
      ],
    );
  }

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
