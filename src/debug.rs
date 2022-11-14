use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use rayon::slice::ParallelSliceMut;

use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use anyhow::Result;

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use crate::util::coloring::ColorMap1d;

#[derive(Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct ColorMap1dRenderer {
  pub width: u32,
  pub height: u32,
  pub filename: String,
  pub color_map: ColorMap1d,
}

impl ColorMap1dRenderer {
  /// Creates a visualization of a color map as a `PNG` image.
  ///
  /// # Errors
  ///
  /// Returns an error, if the generated `PNG` image could not be saved
  /// to disk.
  ///
  pub fn create(&self) -> Result<()> {
    let (w, h) = (self.width as usize, self.height as usize);

    let num_pixel = w * h;

    let mut buf = vec![0_u8; num_pixel * 3];

    let pixel_created = Arc::new(AtomicUsize::new(0));

    buf
      .par_chunks_exact_mut(3)
      .enumerate()
      .for_each(|(i, pixel)| {
        let x = (i % w) as f64;

        let rgb = self.color_map.color(x / w as f64);

        pixel[0] = rgb.r();
        pixel[1] = rgb.g();
        pixel[2] = rgb.b();

        let pc = pixel_created.fetch_add(1, Ordering::SeqCst);

        print!(
          "{}/{} pixels created ({:.2}%)\r",
          pc,
          num_pixel,
          (pc as f32 / num_pixel as f32) * 100.,
        );
      });

    image::save_buffer(
      &self.filename,
      &buf,
      self.width,
      self.height,
      image::ColorType::Rgb8,
    )?;

    println!("\nsuccessfully written: {}", self.filename);

    Ok(())
  }
}
