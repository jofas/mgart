use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use num_complex::Complex64;

use map_macro::vec_no_clone;

use anyhow::Result;

use std::sync::atomic::{AtomicU64, Ordering};

use crate::util::coloring::ColorMap1d;
use crate::util::post_processing::PostProcessing;
use crate::util::sampler::{Sampler, Sampling};
use crate::util::viewport::Viewport;
use crate::util::{print_progress, ComplexNumber};

#[derive(Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct Args {
  pub width: usize,
  pub height: usize,
  pub center: ComplexNumber,
  pub zoom: f64,
  pub iter: u64,
  pub rotation: Option<usize>,
  pub filename: String,
  pub color_map: ColorMap1d,
  pub exponent: f64,
  pub sample_count: u64,
  pub sampler: Sampler,
  #[serde(default)]
  pub post_processing: Vec<PostProcessing>,
}

/// Creates a rendering of a
/// [buddhabrot](https://en.wikipedia.org/wiki/Buddhabrot) as a `PNG`
/// image.
///
/// # Errors
///
/// Returns an error, if the generated `PNG` image could not be saved
/// to disk.
///
pub fn buddhabrot(args: Args) -> Result<()> {
  let aspect_ratio = args.width as f64 / args.height as f64;

  let vp_width = aspect_ratio / args.zoom;
  let vp_height = 1. / args.zoom;

  let grid_delta_x = vp_width / args.width as f64;
  let grid_delta_y = vp_height / args.height as f64;

  let viewport = Viewport::from_center(
    args.center.into(),
    vp_width,
    vp_height,
    grid_delta_x,
    grid_delta_y,
    args.rotation.unwrap_or(0),
  );

  let sampler = args.sampler.distribution(&|c| {
    let (iter, passed_viewport) =
      iter_mandel_check_vp(*c, args.iter, args.exponent, &viewport);

    if iter != args.iter && passed_viewport {
      iter as f64 / args.iter as f64
    } else {
      0.
    }
  });

  println!("starting buddhabrot generation");

  let num_pixel = args.width * args.height;

  let buffer = vec_no_clone![AtomicU64::new(0); num_pixel];

  let processed_samples = AtomicU64::new(0);

  (0..args.sample_count).into_par_iter().for_each(|_| {
    let c = sampler.sample();

    let (j, passed_viewport) =
      iter_mandel_check_vp(c, args.iter, args.exponent, &viewport);

    if j != args.iter && passed_viewport {
      let mut z = c;

      for _ in 0..=j {
        let idx = viewport.grid_pos(&z);

        if let Some((x, y)) = idx {
          buffer[y * args.width + x].fetch_add(1, Ordering::Relaxed);
        }

        z = z.powf(args.exponent) + c;
      }
    }

    let ps = processed_samples.fetch_add(1, Ordering::SeqCst);
    print_progress(ps, args.sample_count, 2500);
  });

  println!("\nbuddhabrot generation finished");

  println!("starting post processing");

  let mut buffer: Vec<f64> =
    buffer.into_iter().map(|x| x.into_inner() as f64).collect();

  for process in args.post_processing {
    process.apply(&mut buffer, args.width, args.height);
  }

  println!("post processing done");

  println!("generating final color values");

  let mut pixels = vec![0_u8; num_pixel * 3];

  pixels
    .chunks_exact_mut(3)
    .enumerate()
    .for_each(|(i, pixel)| {
      let rgb = args.color_map.color(buffer[i]);

      pixel[0] = rgb.r();
      pixel[1] = rgb.g();
      pixel[2] = rgb.b();
    });

  image::save_buffer(
    &args.filename,
    &pixels,
    args.width as u32,
    args.height as u32,
    image::ColorType::Rgb8,
  )?;

  println!("successfully written: {}", args.filename);

  Ok(())
}

fn iter_mandel_check_vp(
  c: Complex64,
  iter: u64,
  exponent: f64,
  viewport: &Viewport,
) -> (u64, bool) {
  let mut z = c;
  let mut z_sqr = z.norm_sqr();

  let mut j = 0;
  let mut passed_viewport = viewport.contains_point(&z);

  while j < iter && z_sqr <= 4.0 {
    z = z.powf(exponent) + c;
    z_sqr = z.norm_sqr();

    if viewport.contains_point(&z) {
      passed_viewport = true;
    }

    j += 1;
  }

  (j, passed_viewport)
}
