use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

pub mod non_local_means;

use non_local_means::NonLocalMeans;

#[derive(
  Serialize, Deserialize, DisplayAsJson, Clone, PartialEq, Debug,
)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "type")]
pub enum Smoothing {
  NonLocalMeans(NonLocalMeans),
}

impl Smoothing {
  pub fn smooth(
    &self,
    buffer: &mut [f64],
    width: usize,
    height: usize,
  ) {
    match self {
      Self::NonLocalMeans(nlm) => nlm.smooth(buffer, width, height),
    }
  }
}
