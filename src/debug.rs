use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use crate::util::coloring::colors::Color;
use crate::util::coloring::ColorMap1d;
use crate::util::frame::Frame;
use crate::util::ProgressPrinter;

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
  #[must_use]
  pub fn create(&self) -> Frame<Color> {
    let (w, h) = (self.width as usize, self.height as usize);

    let num_pixel = w * h;

    let mut frame = Frame::filled_default(w, h);

    let pp = ProgressPrinter::new(num_pixel as u64, 2500);

    frame.par_for_each_mut(|(i, pixel)| {
      let x = (i % w) as f64;

      *pixel = self.color_map.color(x / w as f64).as_color();

      pp.increment();
    });

    frame
  }
}
