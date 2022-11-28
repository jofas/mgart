use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

use crate::util::frame::Frame;

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
  pub fn smooth(&self, frame: &mut Frame<f64>) {
    match self {
      Self::NonLocalMeans(nlm) => nlm.smooth(frame),
    }
  }
}
