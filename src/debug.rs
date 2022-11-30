use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use crate::util::coloring::colors::Color;
use crate::util::coloring::ColorMap1d;
use crate::util::frame::Frame;
use crate::util::ProgressPrinter;

#[derive(Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct ColorMap1dRenderer {
  width: u32,
  height: u32,
  color_map: ColorMap1d,
}

impl ColorMap1dRenderer {
  /// Transforms `self` into a [`Creator`] that can be used
  /// to create an image that shows the color map.
  ///
  #[must_use]
  pub fn creator(self) -> Creator {
    Creator::new(self)
  }
}

pub struct Creator {
  args: ColorMap1dRenderer,
}

impl Creator {
  /// Creates a new instance of [`Creator`].
  ///
  #[must_use]
  pub fn new(args: ColorMap1dRenderer) -> Self {
    Self { args }
  }

  /// Creates a visualization of a color map as a `PNG` image.
  ///
  #[must_use]
  pub fn create(&self) -> Frame<Color> {
    let (w, h) =
      (self.args.width as usize, self.args.height as usize);

    let num_pixel = w * h;

    let mut frame = Frame::filled_default(w, h);

    let pp = ProgressPrinter::new(num_pixel as u64, 2500);

    frame.par_for_each_mut(|(i, pixel)| {
      let x = (i % w) as f64;

      *pixel = self.args.color_map.color(x / w as f64).as_color();

      pp.increment();
    });

    frame
  }
}
