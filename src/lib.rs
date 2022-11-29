#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_lossless)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::float_cmp)]

use serde::Deserialize;

use log::debug;

pub mod buddhabrot;
pub mod debug;
pub mod julia_set;
pub mod util;

use crate::buddhabrot::Buddhabrot;
use crate::debug::ColorMap1dRenderer;
use crate::julia_set::JuliaSet;

#[derive(Deserialize)]
pub struct Algorithm {
  #[serde(flatten)]
  algorithm: AlgorithmInner,
  filename: String,
}

#[derive(Deserialize)]
#[serde(tag = "algorithm")]
#[serde(rename_all = "snake_case")]
pub enum AlgorithmInner {
  JuliaSet(JuliaSet),
  Buddhabrot(Buddhabrot),
  #[serde(rename = "debug.color_map_1d")]
  ColorMap1dRenderer(ColorMap1dRenderer),
}

impl Algorithm {
  /// Executes the rendering process for the given [`Algorithm`],
  /// creating a media file containing the generated artwork.
  ///
  /// # Panics
  ///
  /// Panics if the rendering process fails.
  /// Rendering processes fail, because saving the generated image to
  /// disk was unsuccessful or because the provided
  /// [`configuration`](Self) is faulty.
  ///
  pub fn create(self) {
    match self.algorithm {
      AlgorithmInner::JuliaSet(j) => {
        debug!("generating julia:\n{}", j);
        j.create().save_as_image(&self.filename);
      }
      AlgorithmInner::Buddhabrot(b) => {
        debug!("generating buddhabrot: \n{}", b);
        b.create().save_as_image(&self.filename);
      }
      AlgorithmInner::ColorMap1dRenderer(c) => {
        debug!("generating 1d color map:\n{}", c);
        c.create().save_as_image(&self.filename);
      }
    }
  }
}

#[derive(Deserialize)]
pub struct Algorithms(Vec<Algorithm>);

impl Algorithms {
  /// Executes each [`Algorithm`] successively.
  ///
  /// Multi-threading is implemented inside the rendering process of
  /// each [`Algorithm`].
  ///
  /// # Panics
  ///
  /// If one of the provided algorithms fails.
  ///
  pub fn create(self) {
    for cmd in self.0 {
      cmd.create();
    }
  }
}
