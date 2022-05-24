use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use crate::util::{ColorMap1d, ComplexNumber, Gradient};

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
  #[serde(default = "default_color_map")]
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
  #[serde(default = "default_color_map")]
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

fn default_color_map() -> ColorMap1d {
  ColorMap1d::new(vec![], Gradient::Linear)
}

fn default_julia_set_filename() -> String {
  "julia_set.png".to_owned()
}

fn default_color_map_filename() -> String {
  "color_map.png".to_owned()
}
