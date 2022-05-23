use clap::Args;

use serde::Deserialize;

use crate::util::{ColorMap1D, ComplexNumber};

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

#[derive(Args, Deserialize)]
pub struct JuliaSetArgs {
  #[clap(long, default_value_t = default_width())]
  #[serde(default = "default_width")]
  pub width: u32,
  #[clap(long, default_value_t = default_height())]
  #[serde(default = "default_height")]
  pub height: u32,
  #[clap(long, default_value_t = default_zoom())]
  #[serde(default = "default_zoom")]
  pub zoom: f32,
  #[clap(long, default_value_t = default_zpx())]
  #[serde(default = "default_zpx")]
  pub zpx: f32,
  #[clap(long, default_value_t = default_zpy())]
  #[serde(default = "default_zpy")]
  pub zpy: f32,
  #[clap(long, default_value_t = default_iter())]
  #[serde(default = "default_iter")]
  pub iter: u32,
  #[clap(long, default_value_t = default_filename())]
  #[serde(default = "default_filename")]
  pub filename: String,
  #[clap(long, default_value_t = default_color_map())]
  #[serde(default = "default_color_map")]
  pub color_map: ColorMap1D,
  #[clap(long)]
  #[serde(default)]
  pub c: Option<ComplexNumber>,
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
  ColorMap1D::new(vec![])
}
