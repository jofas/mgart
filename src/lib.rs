#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_lossless)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::float_cmp)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::similar_names)]

use serde::Deserialize;

use anyhow::Result;

use log::debug;

pub mod buddhabrot;
pub mod debug;
pub mod julia_set;
pub mod util;

use crate::buddhabrot::Buddhabrot;
use crate::debug::ColorMap1dRenderer;
use crate::julia_set::JuliaSet;

#[derive(Deserialize)]
#[serde(tag = "algorithm")]
#[serde(rename_all = "snake_case")]
pub enum Algorithm {
  JuliaSet(JuliaSet),
  Buddhabrot(Buddhabrot),
  #[serde(rename = "debug.color_map_1d")]
  ColorMap1dRenderer(ColorMap1dRenderer),
}

impl Algorithm {
  /// Executes the rendering process for the given [`Algorithm`],
  /// creating a media file containing the generated artwork.
  ///
  /// # Errors
  ///
  /// Returns an error if the rendering process fails.
  /// Rendering processes fail, because saving the generated image to
  /// disk was unsuccessful.
  ///
  pub fn create(self) -> Result<()> {
    match self {
      Self::JuliaSet(j) => {
        debug!("generating julia:\n{}", j);
        j.create()
      }
      Self::Buddhabrot(b) => {
        debug!("generating buddhabrot: \n{}", b);
        b.create()
      }
      Self::ColorMap1dRenderer(c) => {
        debug!("generating 1d color map:\n{}", c);
        c.create()
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
  /// # Errors
  ///
  /// If one of the provided algorithms fails, execution is stopped
  /// and the error of the failing alogrithm is returned.
  ///
  pub fn create(self) -> Result<()> {
    for cmd in self.0 {
      cmd.create()?;
    }

    Ok(())
  }
}
