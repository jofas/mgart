use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use num_complex::Complex64;

use map_macro::vec_no_clone;

use log::info;

use std::sync::atomic::{AtomicU64, Ordering};

use crate::util::coloring::colors::Color;
use crate::util::coloring::ColorMap1d;
use crate::util::frame::Frame;
use crate::util::post_processing::PostProcessing;
use crate::util::sampler::{Sampler, Sampling};
use crate::util::viewport::Viewport;
use crate::util::{ComplexNumber, ProgressPrinter};

#[derive(Serialize, Deserialize, DisplayAsJsonPretty, Clone)]
pub struct Buddhabrot {
  width: usize,
  height: usize,
  center: ComplexNumber,
  zoom: f64,
  iter: u64,
  rotation: Option<usize>,
  color_map: ColorMap1d,
  exponent: f64,
  sample_count: u64,
  sampler: Sampler,
  post_processing: Vec<PostProcessing>,
}

impl Buddhabrot {
  /// Creates a new instance of [`Buddhabrot`].
  ///
  #[must_use]
  #[allow(clippy::too_many_arguments)]
  pub fn new(
    width: usize,
    height: usize,
    center: ComplexNumber,
    zoom: f64,
    iter: u64,
    rotation: Option<usize>,
    color_map: ColorMap1d,
    exponent: f64,
    sample_count: u64,
    sampler: Sampler,
    post_processing: Vec<PostProcessing>,
  ) -> Self {
    Self {
      width,
      height,
      center,
      zoom,
      iter,
      rotation,
      color_map,
      exponent,
      sample_count,
      sampler,
      post_processing,
    }
  }

  /// Creates a rendering of a [buddhabrot][buddhabrot] as a `PNG`
  /// image.
  ///
  /// # Errors
  ///
  /// Returns an error, if the configuration is faulty or if the
  /// generated `PNG` image could not be saved to disk.
  ///
  /// [buddhabrot]: https://en.wikipedia.org/wiki/Buddhabrot
  ///
  #[must_use]
  pub fn create(self) -> Frame<Color> {
    let aspect_ratio = self.width as f64 / self.height as f64;

    let vp_width = aspect_ratio / self.zoom;
    let vp_height = 1. / self.zoom;

    let grid_delta_x = vp_width / self.width as f64;
    let grid_delta_y = vp_height / self.height as f64;

    let viewport = Viewport::from_center(
      self.center.into(),
      vp_width,
      vp_height,
      grid_delta_x,
      grid_delta_y,
      self.rotation.unwrap_or(0),
    );

    let sampler = self.sampler.distribution(&|c| {
      let (iter, passed_viewport) = Self::iter_mandel_check_vp(
        *c,
        self.iter,
        self.exponent,
        &viewport,
      );

      if iter != self.iter && passed_viewport {
        iter as f64 / self.iter as f64
      } else {
        0.
      }
    });

    info!("starting buddhabrot generation");

    let buf =
      vec_no_clone![AtomicU64::new(0); self.width * self.height];

    let frame = Frame::new(buf, self.width, self.height);

    let pp = ProgressPrinter::new(self.sample_count, 2500);

    (0..self.sample_count).into_par_iter().for_each(|_| {
      let c = sampler.sample();

      let (j, passed_viewport) = Self::iter_mandel_check_vp(
        c,
        self.iter,
        self.exponent,
        &viewport,
      );

      if j != self.iter && passed_viewport {
        let mut z = c;

        for _ in 0..=j {
          let idx = viewport.rotated_grid_pos(&z);

          if let Some((x, y)) = idx {
            frame[(x, y)].fetch_add(1, Ordering::Relaxed);
          }

          z = z.powf(self.exponent) + c;
        }
      }

      pp.increment();
    });

    info!("buddhabrot generation finished");

    info!("starting post processing");

    let mut frame = frame.map(|x| x.into_inner() as f64);

    for process in self.post_processing {
      process.apply(&mut frame);
    }

    info!("post processing done");

    info!("generating final color values");

    let mut pixels: Frame<Color> =
      Frame::filled_default(self.width, self.height);

    for i in 0..pixels.len() {
      pixels[i] = self.color_map.color(frame[i]).as_color();
    }

    pixels
  }

  fn iter_mandel_check_vp(
    c: Complex64,
    iter: u64,
    exponent: f64,
    viewport: &Viewport,
  ) -> (u64, bool) {
    let mut z = c;
    let mut z_sqr = z.norm_sqr();

    let mut j = 0;
    let mut passed_viewport = viewport.contains_rotated_point(&z);

    while j < iter && z_sqr <= 4.0 {
      z = z.powf(exponent) + c;
      z_sqr = z.norm_sqr();

      if viewport.contains_rotated_point(&z) {
        passed_viewport = true;
      }

      j += 1;
    }

    (j, passed_viewport)
  }
}
