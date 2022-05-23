use clap::Args;

use serde::Deserialize;

use crate::util::{ColorMap1D, RgbaColor};

#[derive(Args, Deserialize)]
pub struct JuliaSetArgs {
  #[clap(default_value = "1920")]
  #[serde(default = "default_width")]
  pub width: u32,
  #[clap(default_value = "1080")]
  #[serde(default = "default_height")]
  pub height: u32,
  #[clap(default_value = "0.5")]
  #[serde(default = "default_zoom")]
  pub zoom: f32,
  #[clap(default_value = "0.")]
  #[serde(default = "default_zpx")]
  pub zpx: f32,
  #[clap(default_value = "0.")]
  #[serde(default = "default_zpy")]
  pub zpy: f32,
  #[clap(default_value = "100")]
  #[serde(default = "default_iter")]
  pub iter: u32,
  #[clap(default_value = "julia_set.png")]
  #[serde(default = "default_filename")]
  pub filename: String,
  #[clap(default_value_t = ColorMap1D::new(vec![]))]
  #[serde(default = "default_color_map")]
  pub color_map: ColorMap1D,
}

fn default_width() -> u32 {
  1920
}

fn default_height() -> u32 {
  1080
}

fn default_zoom() -> f32 {
  0.5
}

fn default_zpx() -> f32 {
  0.
}

fn default_zpy() -> f32 {
  0.
}

fn default_iter() -> u32 {
  100
}

fn default_filename() -> String {
  "julia_set.png".to_owned()
}

fn default_color_map() -> ColorMap1D {
  ColorMap1D::new(vec![
    RgbaColor::new_hex(0xFFFFFFFF),
    RgbaColor::new_hex(0x000000FF),
  ])
}
