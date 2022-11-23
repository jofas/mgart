use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use num_complex::Complex64;

use map_macro::vec_no_clone;

use anyhow::Result;

use log::info;

use std::sync::atomic::{AtomicU64, Ordering};

use crate::util::coloring::ColorMap1d;
use crate::util::post_processing::PostProcessing;
use crate::util::sampler::{Sampler, Sampling};
use crate::util::viewport::Viewport;
use crate::util::{ComplexNumber, ProgressPrinter};

#[derive(Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct Buddhabrot {
  pub width: usize,
  pub height: usize,
  pub center: ComplexNumber,
  pub zoom: f64,
  pub iter: u64,
  pub rotation: Option<usize>,
  pub filename: String,
  pub color_map: ColorMap1d,
  pub exponent: f64,
  pub sample_count: u64,
  pub sampler: Sampler,
  pub post_processing: Vec<PostProcessing>,
}

impl Buddhabrot {
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
  pub fn create(self) -> Result<()> {
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

    let num_pixel = self.width * self.height;

    let buffer = vec_no_clone![AtomicU64::new(0); num_pixel];

    let pp = ProgressPrinter::new(num_pixel as u64, 2500);

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
          let idx = viewport.grid_pos(&z);

          if let Some((x, y)) = idx {
            buffer[y * self.width + x]
              .fetch_add(1, Ordering::Relaxed);
          }

          z = z.powf(self.exponent) + c;
        }
      }

      pp.increment();
    });

    info!("buddhabrot generation finished");

    info!("starting post processing");

    let mut buffer: Vec<f64> =
      buffer.into_iter().map(|x| x.into_inner() as f64).collect();

    for process in self.post_processing {
      process.apply(&mut buffer, self.width, self.height)?;
    }

    info!("post processing done");

    info!("generating final color values");

    let mut pixels = vec![0_u8; num_pixel * 3];

    pixels
      .chunks_exact_mut(3)
      .enumerate()
      .for_each(|(i, pixel)| {
        let rgb = self.color_map.color(buffer[i]);

        pixel[0] = rgb.r();
        pixel[1] = rgb.g();
        pixel[2] = rgb.b();
      });

    image::save_buffer(
      &self.filename,
      &pixels,
      self.width as u32,
      self.height as u32,
      image::ColorType::Rgb8,
    )?;

    info!("successfully written: {}", self.filename);

    Ok(())
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
    let mut passed_viewport = viewport.contains_point(&z);

    while j < iter && z_sqr <= 4.0 {
      z = z.powf(exponent) + c;
      z_sqr = z.norm_sqr();

      if viewport.contains_point(&z) {
        passed_viewport = true;
      }

      j += 1;
    }

    (j, passed_viewport)
  }
}
