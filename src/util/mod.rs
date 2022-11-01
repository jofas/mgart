use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefMutIterator,
  ParallelIterator,
};

use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

use num_complex::Complex64;

use std::f64::consts::PI;
use std::sync::atomic::{AtomicU64, Ordering};

pub mod colors;
pub mod sampler;
pub mod viewport;

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
  Clahe(CLAHE),
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
      Self::Clahe(c) => {
        c.apply(buffer, width, height);
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

        let processed_pixels = AtomicU64::new(0);

        buffer.par_iter_mut().enumerate().for_each(|(i, pixel)| {
          let x = i % width;
          let y = i / width;

          let (wx0, wy0, wx1, wy1) = discrete_bounded_square(
            x,
            y,
            *window_size,
            width,
            height,
          );

          let (x0, y0, x1, y1) =
            discrete_bounded_square(x, y, *n, width, height);

          let c00 = sat[y0 * width + x0];
          let c01 = sat[y0 * width + x1];
          let c10 = sat[y1 * width + x0];
          let c11 = sat[y1 * width + x1];

          let bp = c00 + c11 - c01 - c10;
          let bp = bp / ((x1 - x0 + 1) * (y1 - y0 + 1)) as f64;

          let mut s = 0.;
          let mut cp = 0.;

          for x in wx0..=wx1 {
            for y in wy0..=wy1 {
              let (x0, y0, x1, y1) =
                discrete_bounded_square(x, y, *n, width, height);

              let c00 = sat[y0 * width + x0];
              let c01 = sat[y0 * width + x1];
              let c10 = sat[y1 * width + x0];
              let c11 = sat[y1 * width + x1];

              let bq = c00 + c11 - c01 - c10;
              let bq = bq / ((x1 - x0 + 1) * (y1 - y0 + 1)) as f64;

              let fpq = (-(bq - bp).powi(2) / h.powi(2)).exp();

              s += *pixel * fpq;
              cp += fpq;
            }
          }

          *pixel = s / cp;

          let pc = processed_pixels.fetch_add(1, Ordering::SeqCst);
          print_progress(pc, num_pixel as u64, 100);
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

/*
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "process")]
pub struct Sampler {
  Uniform,
  Kde,
}

impl Sampler {
  pub fn into_sampler(self) -> dyn Sampling {
    match Self {
      Kde => {
       samples =
      }
    }
  }
}
*/

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct CLAHE {
  contrast_limit: usize,
  bin_count: usize,
  tile_size_x: usize,
  tile_size_y: usize,
}

impl CLAHE {
  pub fn new(
    contrast_limit: usize,
    bin_count: usize,
    tile_size_x: usize,
    tile_size_y: usize,
  ) -> Self {
    Self {
      contrast_limit,
      bin_count,
      tile_size_x,
      tile_size_y,
    }
  }

  pub fn apply(
    &self,
    buffer: &mut [f64],
    width: usize,
    height: usize,
  ) {
    if width % self.tile_size_x != 0 || height % self.tile_size_y != 0
    {
      panic!("width and height must be divisible by tile_size_x and tile_size_y, respectively");
    }

    let tiles_w = width / self.tile_size_x;
    let tiles_h = height / self.tile_size_y;

    let tiles = self.tiles(buffer, width, tiles_w, tiles_h);

    for (i, v) in buffer.iter_mut().enumerate() {
      let x = i % width;
      let y = i / width;

      let pos = Pos::new(
        x % self.tile_size_x,
        y % self.tile_size_y,
        self.tile_size_x,
        self.tile_size_y,
      );

      let (x_tile, y_tile) = self.tile_indices(x, y);

      if let Pos::Center = pos {
        let i = y_tile as usize * tiles_w + x_tile as usize;
        *v = tiles[i].transform(*v);
        continue;
      }

      let (nw, ne, se, sw) = self.interpolation_tiles(
        &pos, x_tile, y_tile, &tiles, tiles_w, tiles_h,
      );

      let (dn, ds, dw, de) = self.interpolation_distances(x, y, &pos);

      let (dn, ds, dw, de) = self.handle_corners_and_borders(
        &nw, &ne, &se, &sw, dn, ds, dw, de,
      );

      let q_nw = nw.transform(*v) * dn * dw;
      let q_ne = ne.transform(*v) * dn * de;
      let q_se = se.transform(*v) * ds * de;
      let q_sw = sw.transform(*v) * ds * dw;

      *v = q_nw + q_ne + q_se + q_sw;
    }
  }

  fn tiles(
    &self,
    buffer: &[f64],
    width: usize,
    tiles_w: usize,
    tiles_h: usize,
  ) -> Vec<Tile> {
    let mut tiles: Vec<Tile> = Vec::with_capacity(tiles_w * tiles_h);

    for start_block in 0..tiles_h {
      for offset in 0..tiles_w {
        let stride = Strided::new(
          width,
          self.tile_size_x,
          offset * self.tile_size_x,
          Some(self.tile_size_y),
          Some(start_block * self.tile_size_y),
          &buffer,
        );

        tiles.push(Tile::new(
          stride,
          self.bin_count,
          self.contrast_limit,
        ));
      }
    }

    tiles
  }

  fn interpolation_tiles<'a>(
    &self,
    pos: &Pos,
    x_tile: isize,
    y_tile: isize,
    tiles: &'a [Tile],
    tiles_w: usize,
    tiles_h: usize,
  ) -> (
    Option<&'a Tile>,
    Option<&'a Tile>,
    Option<&'a Tile>,
    Option<&'a Tile>,
  ) {
    let (x1, x2, y1, y2) = match pos {
      Pos::NW => (x_tile - 1, x_tile, y_tile - 1, y_tile),
      Pos::NE => (x_tile, x_tile + 1, y_tile - 1, y_tile),
      Pos::SE => (x_tile, x_tile + 1, y_tile, y_tile + 1),
      Pos::SW => (x_tile - 1, x_tile, y_tile, y_tile + 1),
      Pos::Center => panic!("Cannot handle pixels at tile center"),
    };

    let mut nw = None;
    let mut ne = None;
    let mut se = None;
    let mut sw = None;

    if x1 >= 0 && y1 >= 0 {
      nw = Some(&tiles[y1 as usize * tiles_w + x1 as usize]);
    }

    if x2 < tiles_w as isize && y1 >= 0 {
      ne = Some(&tiles[y1 as usize * tiles_w + x2 as usize]);
    }

    if x2 < tiles_w as isize && y2 < tiles_h as isize {
      se = Some(&tiles[y2 as usize * tiles_w + x2 as usize]);
    }

    if x1 >= 0 && y2 < tiles_h as isize {
      sw = Some(&tiles[y2 as usize * tiles_w + x1 as usize]);
    }

    (nw, ne, se, sw)
  }

  fn interpolation_distances(
    &self,
    x: usize,
    y: usize,
    pos: &Pos,
  ) -> (f64, f64, f64, f64) {
    let center_x = self.tile_size_x as f64 / 2.;
    let center_y = self.tile_size_y as f64 / 2.;

    let dx = (center_x - (x % self.tile_size_x) as f64).abs();
    let dy = (center_y - (y % self.tile_size_y) as f64).abs();

    let dx = dx / self.tile_size_x as f64;
    let dy = dy / self.tile_size_y as f64;

    match pos {
      Pos::NW => (dy, 1. - dy, dx, 1. - dx),
      Pos::NE => (dy, 1. - dy, 1. - dx, dx),
      Pos::SE => (1. - dy, dy, 1. - dx, dx),
      Pos::SW => (1. - dy, dy, dx, 1. - dx),
      Pos::Center => panic!("Unable to handle pixels at tile center"),
    }
  }

  fn handle_corners_and_borders(
    &self,
    nw: &Option<&Tile>,
    ne: &Option<&Tile>,
    se: &Option<&Tile>,
    sw: &Option<&Tile>,
    dn: f64,
    ds: f64,
    dw: f64,
    de: f64,
  ) -> (f64, f64, f64, f64) {
    match (nw, ne, se, sw) {
      // corners
      (Some(_), None, None, None) => (1., 0., 1., 0.),
      (None, Some(_), None, None) => (1., 0., 0., 1.),
      (None, None, Some(_), None) => (0., 1., 0., 1.),
      (None, None, None, Some(_)) => (0., 1., 1., 0.),
      // borders
      (Some(_), Some(_), None, None) => (1., 0., dw, de),
      (Some(_), None, None, Some(_)) => (dn, ds, 1., 0.),
      (None, Some(_), Some(_), None) => (dn, ds, 1., 0.),
      (None, None, Some(_), Some(_)) => (1., 0., dw, de),
      // center
      (Some(_), Some(_), Some(_), Some(_)) => (dn, ds, dw, de),
      _ => panic!("impossible state"),
    }
  }

  fn tile_indices(&self, x: usize, y: usize) -> (isize, isize) {
    (
      (x / self.tile_size_x) as isize,
      (y / self.tile_size_y) as isize,
    )
  }
}

#[derive(Debug)]
enum Pos {
  NW,
  NE,
  SE,
  SW,
  Center,
}

impl Pos {
  fn new(x: usize, y: usize, x_max: usize, y_max: usize) -> Self {
    let pos_x = PosH::new(x, x_max);
    let pos_y = PosV::new(y, y_max);

    match (pos_y, pos_x) {
      (PosV::N, PosH::W) => Self::NW,
      (PosV::N, PosH::E) => Self::NE,
      (PosV::N, PosH::Center) => Self::NW,
      (PosV::S, PosH::E) => Self::SE,
      (PosV::S, PosH::W) => Self::SW,
      (PosV::S, PosH::Center) => Self::SE,
      (PosV::Center, PosH::E) => Self::NE,
      (PosV::Center, PosH::W) => Self::SW,
      (PosV::Center, PosH::Center) => Self::Center,
    }
  }
}

enum PosH {
  W,
  E,
  Center,
}

impl PosH {
  fn new(x: usize, x_max: usize) -> Self {
    if x_max % 2 == 0 {
      if x < x_max / 2 {
        Self::W
      } else {
        Self::E
      }
    } else {
      if x < x_max / 2 {
        Self::W
      } else if x > x_max / 2 {
        Self::E
      } else {
        Self::Center
      }
    }
  }
}

enum PosV {
  N,
  S,
  Center,
}

impl PosV {
  fn new(y: usize, y_max: usize) -> Self {
    if y_max % 2 == 0 {
      if y < y_max / 2 {
        Self::N
      } else {
        Self::S
      }
    } else {
      if y < y_max / 2 {
        Self::N
      } else if y > y_max / 2 {
        Self::S
      } else {
        Self::Center
      }
    }
  }
}

trait HistogramEqualization {
  fn transform(&self, v: f64) -> f64;
}

struct Tile {
  hist: Vec<usize>,
  cdf_min: usize,
  n: usize,
}

impl Tile {
  pub fn new(
    buffer: impl Iterator<Item = f64>,
    bin_count: usize,
    contrast_limit: usize,
  ) -> Self {
    let mut hist = vec![0; bin_count];
    let mut n = 0;

    // contrast limiting
    let mut clv = 0;

    for v in buffer {
      let bin = (v * (bin_count - 1) as f64) as usize;

      if hist[bin] < contrast_limit {
        hist[bin] += 1;
      } else {
        clv += 1;
      }

      n += 1;
    }

    loop {
      let clv_bin = clv / bin_count;

      let mut new_clv = 0;

      for b in &mut hist {
        *b += clv_bin;

        if *b > contrast_limit {
          new_clv += *b - contrast_limit;
          *b = contrast_limit;
        }
      }

      if new_clv == 0 || new_clv == clv {
        break;
      }

      clv = new_clv;
    }

    // sum up histograms
    for i in 1..hist.len() {
      hist[i] += hist[i - 1];
    }

    // minimum value in histogram
    let mut cdf_min = 0;
    for b in &hist {
      if *b > 0 {
        cdf_min = *b;
        break;
      }
    }

    Self { hist, cdf_min, n }
  }
}

impl HistogramEqualization for Tile {
  fn transform(&self, v: f64) -> f64 {
    let bin = (v * (self.hist.len() - 1) as f64) as usize;

    let res = (self.hist[bin] - self.cdf_min) as f64;

    res / (self.n - self.cdf_min) as f64
  }
}

impl HistogramEqualization for Option<&Tile> {
  fn transform(&self, v: f64) -> f64 {
    match self {
      Some(t) => t.transform(v),
      None => 0.,
    }
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

impl<'a, T: Copy> Iterator for Strided<'a, T> {
  type Item = T;

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

    self.buffer.get(i).map(|x| *x)
  }
}

pub fn discrete_bounded_square(
  x: usize,
  y: usize,
  n: usize,
  width: usize,
  height: usize,
) -> (usize, usize, usize, usize) {
  let nh = n / 2;

  let xupper = (x + nh).min(width - 1);
  let yupper = (y + nh).min(height - 1);

  let (x, y) = if n % 2 == 0 { (x + 1, y + 1) } else { (x, y) };

  let xlower = x.checked_sub(nh).unwrap_or(0);
  let ylower = y.checked_sub(nh).unwrap_or(0);

  (xlower, ylower, xupper, yupper)
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

pub fn print_progress(i: u64, n: u64, interval: u64) {
  if i % interval == interval - 1 || i == n - 1 {
    let p = i as f64 / n as f64 * 100.;
    print!("{}/{} iterations done ({:.2}%)\r", i + 1, n, p);
  }
}

#[cfg(test)]
mod tests {
  use super::{
    discrete_bounded_square, summed_area_table, Gradient, Strided,
    CLAHE,
  };

  #[test]
  fn tile_indexing() {
    let c = CLAHE::new(0, 0, 8, 8);

    assert_eq!((0, 0), c.tile_indices(0, 0));
    assert_eq!((0, 0), c.tile_indices(7, 0));
    assert_eq!((0, 0), c.tile_indices(0, 7));
    assert_eq!((0, 0), c.tile_indices(7, 7));

    assert_eq!((1, 1), c.tile_indices(8, 8));
    assert_eq!((1, 1), c.tile_indices(15, 8));
    assert_eq!((1, 1), c.tile_indices(8, 15));
    assert_eq!((1, 1), c.tile_indices(15, 15));

    assert_eq!((2, 2), c.tile_indices(16, 16));
    assert_eq!((2, 2), c.tile_indices(23, 16));
    assert_eq!((2, 2), c.tile_indices(16, 23));
    assert_eq!((2, 2), c.tile_indices(23, 23));

    assert_eq!((3, 3), c.tile_indices(24, 24));
    assert_eq!((3, 3), c.tile_indices(31, 24));
    assert_eq!((3, 3), c.tile_indices(24, 31));
    assert_eq!((3, 3), c.tile_indices(31, 31));
  }

  #[test]
  fn strided() {
    let buf: Vec<f64> = (0..16).map(|x| x as f64).collect();

    let s = Strided::new(4, 2, 0, None, None, &buf);
    let v: Vec<f64> = s.map(|x| x).collect();
    assert_eq!(v, vec![0., 1., 4., 5., 8., 9., 12., 13.]);

    let s = Strided::new(4, 2, 2, None, None, &buf);
    let v: Vec<f64> = s.map(|x| x).collect();
    assert_eq!(v, vec![2., 3., 6., 7., 10., 11., 14., 15.]);

    let s = Strided::new(4, 3, 1, None, None, &buf);
    let v: Vec<f64> = s.map(|x| x).collect();
    assert_eq!(
      v,
      vec![1., 2., 3., 5., 6., 7., 9., 10., 11., 13., 14., 15.]
    );

    let s = Strided::new(4, 2, 0, Some(2), None, &buf);
    let v: Vec<f64> = s.map(|x| x).collect();
    assert_eq!(v, vec![0., 1., 4., 5.]);

    let s = Strided::new(4, 2, 2, Some(2), None, &buf);
    let v: Vec<f64> = s.map(|x| x).collect();
    assert_eq!(v, vec![2., 3., 6., 7.]);

    let s = Strided::new(4, 2, 0, Some(2), Some(2), &buf);
    let v: Vec<f64> = s.map(|x| x).collect();
    assert_eq!(v, vec![8., 9., 12., 13.]);

    let s = Strided::new(4, 2, 2, Some(2), Some(2), &buf);
    let v: Vec<f64> = s.map(|x| x).collect();
    assert_eq!(v, vec![10., 11., 14., 15.]);
  }

  #[test]
  fn dbs() {
    let res = discrete_bounded_square(0, 0, 3, 4, 4);
    assert_eq!(res, (0, 0, 1, 1));

    let res = discrete_bounded_square(1, 1, 3, 4, 4);
    assert_eq!(res, (0, 0, 2, 2));

    let res = discrete_bounded_square(2, 2, 3, 4, 4);
    assert_eq!(res, (1, 1, 3, 3));

    let res = discrete_bounded_square(3, 3, 3, 4, 4);
    assert_eq!(res, (2, 2, 3, 3));

    let res = discrete_bounded_square(0, 0, 2, 4, 4);
    assert_eq!(res, (0, 0, 1, 1));

    let res = discrete_bounded_square(1, 1, 2, 4, 4);
    assert_eq!(res, (1, 1, 2, 2));

    let res = discrete_bounded_square(2, 2, 2, 4, 4);
    assert_eq!(res, (2, 2, 3, 3));

    let res = discrete_bounded_square(3, 3, 2, 4, 4);
    assert_eq!(res, (3, 3, 3, 3));
  }

  #[test]
  fn sat() {
    let image = [1.; 16];
    let width = 4;

    let sat = summed_area_table(&image, width);

    assert_eq!(
      sat,
      vec![
        [1., 2., 3., 4.],
        [2., 4., 6., 8.],
        [3., 6., 9., 12.],
        [4., 8., 12., 16.]
      ]
      .into_iter()
      .flatten()
      .collect::<Vec<f64>>(),
    );

    let image: Vec<f64> = (1..=16).map(|x| x as f64).collect();
    let width = 4;

    let sat = summed_area_table(&image, width);

    assert_eq!(
      sat,
      vec![
        [1., 3., 6., 10.],
        [6., 14., 24., 36.],
        [15., 33., 54., 78.],
        [28., 60., 96., 136.]
      ]
      .into_iter()
      .flatten()
      .collect::<Vec<f64>>(),
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
}
