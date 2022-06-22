use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use crate::util::{ColorMap1d, ComplexNumber};

#[derive(Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct BuddhabrotArgs {
  #[serde(default = "default_width")]
  pub width: usize,
  #[serde(default = "default_height")]
  pub height: usize,
  #[serde(default = "default_zoom")]
  pub zoom: f64,
  #[serde(default = "default_zpx")]
  pub zpx: f64,
  #[serde(default = "default_zpy")]
  pub zpy: f64,
  #[serde(default = "default_iter")]
  pub iter: u32,
  #[serde(default = "default_buddhabrot_filename")]
  pub filename: String,
  #[serde(default)]
  pub color_map: ColorMap1d,
  #[serde(default = "default_sample_count")]
  pub sample_count: u64,
  #[serde(default = "default_deviation_factor")]
  pub deviation_factor: f64,
  #[serde(default = "default_transition_threshold")]
  pub transition_threshold: f64,
}

#[derive(Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct JuliaSetArgs {
  #[serde(default = "default_width")]
  pub width: usize,
  #[serde(default = "default_height")]
  pub height: usize,
  #[serde(default = "default_zoom")]
  pub zoom: f64,
  #[serde(default = "default_zpx")]
  pub zpx: f64,
  #[serde(default = "default_zpy")]
  pub zpy: f64,
  #[serde(default = "default_iter")]
  pub iter: u32,
  #[serde(default = "default_julia_set_filename")]
  pub filename: String,
  #[serde(default)]
  pub color_map: ColorMap1d,
  #[serde(default)]
  pub c: Option<ComplexNumber>,
}

#[derive(Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct ColorMap1dArgs {
  #[serde(default = "default_width")]
  pub width: usize,
  #[serde(default = "default_height")]
  pub height: usize,
  #[serde(default = "default_color_map_filename")]
  pub filename: String,
  #[serde(default)]
  pub color_map: ColorMap1d,
}

fn default_width() -> usize {
  1920
}

fn default_height() -> usize {
  1080
}

fn default_zoom() -> f64 {
  0.5
}

fn default_zpx() -> f64 {
  0.
}

fn default_zpy() -> f64 {
  0.
}

fn default_iter() -> u32 {
  100
}

fn default_sample_count() -> u64 {
  100_000_000
}

fn default_deviation_factor() -> f64 {
  0.1
}

fn default_transition_threshold() -> f64 {
  0.1
}

fn default_buddhabrot_filename() -> String {
  "buddhabrot.png".to_owned()
}

fn default_julia_set_filename() -> String {
  "julia_set.png".to_owned()
}

fn default_color_map_filename() -> String {
  "color_map.png".to_owned()
}
