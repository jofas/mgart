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
  /// Creates a new instance of [`ColorMap1d`].
  ///
  #[must_use]
  pub fn new(width: u32, height: u32, color_map: ColorMap1d) -> Self {
    Self {
      width,
      height,
      color_map,
    }
  }

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
  /// # Panics
  ///
  /// Panics, if the configuration is faulty, e.g if there happens an
  /// overflow.
  ///
  #[must_use]
  pub fn create(&self) -> Frame<Color> {
    let mut frame =
      Frame::filled_default(self.args.width, self.args.height);

    let pp = ProgressPrinter::new(frame.len() as u64, 2500);

    frame.par_for_each_mut(|(i, pixel)| {
      let x = f64::from(u32::try_from(i).unwrap() % self.args.width);

      *pixel = self
        .args
        .color_map
        .color(x / f64::from(self.args.width))
        .as_color();

      pp.increment();
    });

    frame
  }
}
