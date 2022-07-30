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

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "process")]
pub enum PostProcessing {
  Normalize,
  Clamp { min: f64, max: f64 },
  ClampAndNormalize { min: f64, max: f64 },
  Gradient(Gradient),
  Smoothing(Smoothing),
}

impl PostProcessing {
  pub fn apply(
    &self,
    buffer: &mut [f64],
    width: usize,
    height: usize,
  ) {
    match self {
      Self::Normalize => {
        let (min, max) = min_max(buffer);

        for v in buffer {
          *v = (*v - min) / (max - min);
        }
      }
      Self::Clamp { min, max } => {
        for v in buffer {
          *v = v.clamp(*min, *max);
        }
      }
      Self::ClampAndNormalize { min, max } => {
        for v in buffer {
          *v = (v.clamp(*min, *max) - min) / (max - min);
        }
      }
      Self::Gradient(g) => {
        for v in buffer {
          *v = g.apply(*v);
        }
      }
      Self::Smoothing(s) => {
        s.smooth(buffer, width, height);
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
        let n = *n as f64;
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
    let f = self.gradient.apply(f);

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

#[derive(
  Serialize, Deserialize, DisplayAsJson, Clone, PartialEq, Debug,
)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "type")]
pub enum Smoothing {
  NonLocalMeans {
    n: usize,
    window_size: usize,
    h: f64,
  },
}

impl Smoothing {
  pub fn smooth(
    &self,
    buffer: &mut [f64],
    width: usize,
    height: usize,
  ) {
    match self {
      Self::NonLocalMeans { n, window_size, h } => {
        let sat = summed_area_table(&buffer, width);

        let num_pixel = buffer.len();

        let processed_pixels = AtomicU32::new(0);

        buffer.par_iter_mut().enumerate().for_each(|(i, pixel)| {
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

              let fpq = (-((bq - bp).powi(2) / h)).exp();

              s += *pixel * fpq;
              cp += fpq;
            }
          }

          *pixel = s / cp;

          let pc = processed_pixels.fetch_add(1, Ordering::SeqCst);
          print_progress(pc, num_pixel as u32, 100);
        });
      }
    }
  }
}

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

// TODO: implement the A (adaptive) part
pub struct CLAHE {
  contrast_limit: usize,
  bin_count: usize,
  tile_size: usize, // TODO: x and y
}

impl CLAHE {
  pub fn new(
    contrast_limit: usize,
    bin_count: usize,
    tile_size: usize,
  ) -> Self {
    Self {
      contrast_limit,
      bin_count,
      tile_size,
    }
  }

  pub fn apply(
    &self,
    buffer: &mut [f64],
    width: usize,
    height: usize,
  ) {
    if width % self.tile_size != 0 || height % self.tile_size != 0 {
      panic!("width and height must be divisible by tile_size");
    }

    // TODO: create histogram for each tile
    //
    // Tile struct containing histogram
    //
    // 3 possilbe pixels:
    //  1. center tile || corner (1 TF)
    //  2. border (2 TF)
    //  3. area (4 TF)

    let bins = self.create_bins(buffer);

    let mut cdf_min = 0;
    for b in &bins {
      if *b > 0 {
        cdf_min = *b;
        break;
      }
    }

    let n = buffer.len();

    for v in buffer {
      let bin = (*v * (self.bin_count - 1) as f64) as usize;

      // transformation function
      *v = (bins[bin] - cdf_min) as f64 / (n - cdf_min) as f64;
    }
  }

  pub(self) fn create_bins(&self, buffer: &[f64]) -> Vec<usize> {
    let mut bins = vec![0; self.bin_count];

    let mut clv = 0;

    for v in buffer.iter() {
      let bin = (*v * (self.bin_count - 1) as f64) as usize;

      if bins[bin] < self.contrast_limit {
        bins[bin] += 1;
      } else {
        clv += 1;
      }
    }

    clv = clv / self.bin_count;

    if clv > 0 {
      for b in &mut bins {
        *b += clv;
      }
    }

    for i in 1..bins.len() {
      bins[i] += bins[i - 1];
    }

    bins
  }
}

struct Tile {}

impl Tile {
  pub fn new(
    buffer: &[f64],
    block_size: usize,
    elements: usize,
    offset: usize,
  ) -> Self {
    assert_eq!(buffer.len() % block_size, 0);
    Self {}
  }
}

struct Strided<'a, T> {
  block_size: usize,
  elements: usize,
  offset: usize,
  block_count: Option<usize>,
  start_block: usize,
  buffer: &'a [T],
  block: usize,
  element: usize,
}

impl<'a, T> Strided<'a, T> {
  pub fn new(
    block_size: usize,
    elements: usize,
    offset: usize,
    block_count: Option<usize>,
    start_block: Option<usize>,
    buffer: &'a [T],
  ) -> Self {
    Self {
      block_size,
      elements,
      offset,
      block_count,
      start_block: start_block.unwrap_or(0),
      buffer,
      block: 0,
      element: 0,
    }
  }
}

impl<'a, T> Iterator for Strided<'a, T> {
  type Item = &'a T;

  fn next(&mut self) -> Option<Self::Item> {
    // early exit if block_count is defined and has been reached
    match self.block_count {
      Some(bc) if bc == self.block => return None,
      _ => (),
    }

    let i = (self.start_block + self.block) * self.block_size
      + self.offset
      + self.element;

    if self.element < self.elements - 1 {
      self.element += 1;
    } else {
      self.element = 0;
      self.block += 1;
    }

    self.buffer.get(i)
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
  use super::{
    grid_pos, summed_area_table, Gradient, Strided, Viewport, CLAHE,
  };

  /*
  #[test]
  fn clahe_bins() {
    let image = vec![0., 0.34, 0.34, 0., 0.67, 1.];

    let c = CLAHE::new(2, 4);

    let bins = c.create_bins(&image);

    assert_eq!(bins, vec![2, 4, 5, 6]);
  }

  #[test]
  fn clahe_cl_no_clip() {
    let image = vec![0., 0.34, 0.34, 0., 0.67, 1.];

    let c = CLAHE::new(1, 4);

    let bins = c.create_bins(&image);

    assert_eq!(bins, vec![1, 2, 3, 4]);
  }

  #[test]
  fn clahe_cl_clip() {
    let image = vec![0., 0.34, 0.34, 0., 0.67, 1., 0., 0.34];

    let c = CLAHE::new(1, 4);

    let bins = c.create_bins(&image);

    assert_eq!(bins, vec![2, 4, 6, 8]);
  }

  #[test]
  fn clahe() {
    let mut image = vec![0., 0.34, 0.34, 0., 0.67, 1.];

    let c = CLAHE::new(2, 4);

    c.apply(&mut image);

    assert_eq!(image, vec![0., 0.5, 0.5, 0., 0.75, 1.]);
  }
  */

  #[test]
  fn strided() {
    let buf: Vec<f64> = (0..16).map(|x| x as f64).collect();

    let s = Strided::new(4, 2, 0, None, None, &buf);
    let v: Vec<f64> = s.map(|x| *x).collect();
    assert_eq!(v, vec![0., 1., 4., 5., 8., 9., 12., 13.]);

    let s = Strided::new(4, 2, 2, None, None, &buf);
    let v: Vec<f64> = s.map(|x| *x).collect();
    assert_eq!(v, vec![2., 3., 6., 7., 10., 11., 14., 15.]);

    let s = Strided::new(4, 3, 1, None, None, &buf);
    let v: Vec<f64> = s.map(|x| *x).collect();
    assert_eq!(
      v,
      vec![1., 2., 3., 5., 6., 7., 9., 10., 11., 13., 14., 15.]
    );

    let s = Strided::new(4, 2, 0, Some(2), None, &buf);
    let v: Vec<f64> = s.map(|x| *x).collect();
    assert_eq!(v, vec![0., 1., 4., 5.]);

    let s = Strided::new(4, 2, 2, Some(2), None, &buf);
    let v: Vec<f64> = s.map(|x| *x).collect();
    assert_eq!(v, vec![2., 3., 6., 7.]);

    let s = Strided::new(4, 2, 0, Some(2), Some(2), &buf);
    let v: Vec<f64> = s.map(|x| *x).collect();
    assert_eq!(v, vec![8., 9., 12., 13.]);

    let s = Strided::new(4, 2, 2, Some(2), Some(2), &buf);
    let v: Vec<f64> = s.map(|x| *x).collect();
    assert_eq!(v, vec![10., 11., 14., 15.]);
  }

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

    assert_eq!(g.apply(0.), 0.);
    assert_eq!(g.apply(0.5), 0.5);
    assert_eq!(g.apply(1.), 1.);

    let g = Gradient::Linear { factor: 2. };

    assert_eq!(g.apply(0.), 0.);
    assert_eq!(g.apply(0.25), 0.5);
    assert_eq!(g.apply(0.5), 1.);
    assert_eq!(g.apply(0.75), 0.5);
    assert_eq!(g.apply(1.), 1.);
  }

  #[test]
  fn discrete_gradient() {
    let g = Gradient::Discrete {
      n: 2,
      gradient: Box::new(Gradient::Linear { factor: 1. }),
    };

    assert_eq!(g.apply(0.), 0.);
    assert_eq!(g.apply(0.25), 0.);
    assert_eq!(g.apply(0.5), 1.);
    assert_eq!(g.apply(0.75), 1.);
    assert_eq!(g.apply(1.), 1.);

    let g = Gradient::Discrete {
      n: 3,
      gradient: Box::new(Gradient::Linear { factor: 1. }),
    };

    assert_eq!(g.apply(0.), 0.);
    assert_eq!(g.apply(0.33), 0.);
    assert_eq!(g.apply(0.34), 0.5);
    assert_eq!(g.apply(0.66), 0.5);
    assert_eq!(g.apply(0.67), 1.);
    assert_eq!(g.apply(1.), 1.);

    let g = Gradient::Discrete {
      n: 4,
      gradient: Box::new(Gradient::Linear { factor: 1. }),
    };

    assert_eq!(g.apply(0.), 0.);
    assert_eq!(g.apply(0.24), 0.);
    assert_eq!((g.apply(0.25) * 100.).floor(), 33.);
    assert_eq!((g.apply(0.49) * 100.).floor(), 33.);
    assert_eq!((g.apply(0.5) * 100.).floor(), 66.);
    assert_eq!((g.apply(0.74) * 100.).floor(), 66.);
    assert_eq!(g.apply(0.75), 1.);
    assert_eq!(g.apply(1.), 1.);
  }

  #[test]
  fn grid_pos1() {
    let vp = Viewport::new(0., 0., 2., 2.);

    let (x, y) = grid_pos(1.5, 0.5, 1., 1., &vp).unwrap();

    assert_eq!(x, 1);
    assert_eq!(y, 0);
  }

  #[test]
  fn grid_pos2() {
    let vp = Viewport::new(0., 0., 2., 2.);

    let p = grid_pos(1.5, -0.5, 1., 1., &vp);

    assert_eq!(p, None);
  }

  #[test]
  fn grid_pos3() {
    let vp = Viewport::new(0., 0., 2., 2.);

    let p = grid_pos(1.5, 2.5, 1., 1., &vp);

    assert_eq!(p, None);
  }

  #[test]
  fn grid_pos4() {
    let vp = Viewport::new(-50., -50., 50., 50.);

    let (x, y) = grid_pos(0., 0., 1., 1., &vp).unwrap();

    assert_eq!(x, 50);
    assert_eq!(y, 50);
  }

  #[test]
  fn grid_pos5() {
    let vp = Viewport::new(-1., -1., 1., 1.);

    let (x, y) = grid_pos(0., 0., 0.01, 0.005, &vp).unwrap();

    assert_eq!(x, 100);
    assert_eq!(y, 200);
  }

  #[test]
  fn grid_pos6() {
    let vp = Viewport::new(-1., -1., 1., 1.);

    let (x, y) = grid_pos(0.999, 0.999, 0.01, 0.005, &vp).unwrap();

    assert_eq!(x, 199);
    assert_eq!(y, 399);
  }
}
