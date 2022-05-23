use clap::Args;

use serde::Deserialize;

use crate::util::ColorMap1D;

#[derive(Args, Deserialize)]
pub struct JuliaSetArgs {
  width: u32,
  height: u32,
  zoom: f32,
  zpx: f32,
  zpy: f32,
  iter: u32,
  filename: String,
  color_map: ColorMap1D,
}
