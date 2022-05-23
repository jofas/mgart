use clap::Args;

use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use crate::util::{ColorMap1d, ColorMethod, ComplexNumber};

#[derive(Deserialize)]
pub struct Config(Vec<JuliaSetArgs>);

impl Config {
  pub fn inner(&self) -> &[JuliaSetArgs] {
    &self.0
  }

  pub fn into_inner(self) -> Vec<JuliaSetArgs> {
    self.0
  }
}

#[derive(Args, Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct JuliaSetArgs {
  #[clap(long, default_value_t = default_width())]
  #[serde(default = "default_width")]
  pub width: usize,
  #[clap(long, default_value_t = default_height())]
  #[serde(default = "default_height")]
  pub height: usize,
  #[clap(long, default_value_t = default_zoom())]
  #[serde(default = "default_zoom")]
  pub zoom: f64,
  #[clap(long, default_value_t = default_zpx())]
  #[serde(default = "default_zpx")]
  pub zpx: f64,
  #[clap(long, default_value_t = default_zpy())]
  #[serde(default = "default_zpy")]
  pub zpy: f64,
  #[clap(long, default_value_t = default_iter())]
  #[serde(default = "default_iter")]
  pub iter: u32,
  #[clap(long, default_value_t = default_julia_set_filename())]
  #[serde(default = "default_julia_set_filename")]
  pub filename: String,
  #[clap(long, default_value_t = default_color_map())]
  #[serde(default = "default_color_map")]
  pub color_map: ColorMap1d,
  #[clap(long, default_value_t = default_color_method())]
  #[serde(default = "default_color_method")]
  pub color_method: ColorMethod,
  #[clap(long)]
  #[serde(default)]
  pub c: Option<ComplexNumber>,
}

#[derive(Args, Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct ColorMap1dArgs {
  #[clap(long, default_value_t = default_width())]
  #[serde(default = "default_width")]
  pub width: usize,
  #[clap(long, default_value_t = default_height())]
  #[serde(default = "default_height")]
  pub height: usize,
  #[clap(long, default_value_t = default_color_map_filename())]
  #[serde(default = "default_color_map_filename")]
  pub filename: String,
  #[clap(long, default_value_t = default_color_map())]
  #[serde(default = "default_color_map")]
  pub color_map: ColorMap1d,
  #[clap(long, default_value_t = default_color_method())]
  #[serde(default = "default_color_method")]
  pub color_method: ColorMethod,
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
  ColorMap1d::new(vec![])
}

fn default_color_method() -> ColorMethod {
  ColorMethod::Linear
}

fn default_julia_set_filename() -> String {
  "julia_set.png".to_owned()
}

fn default_color_map_filename() -> String {
  "color_map.png".to_owned()
}
