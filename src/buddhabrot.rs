use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use num_complex::Complex64;

use rand::random;

use map_macro::vec_no_clone;

use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};

use crate::util::sampler::{Sampler as Sampler_, Uniform, KDE};
use crate::util::{
  grid_pos, print_progress, ColorMap1d, PostProcessing, Viewport,
};

#[derive(Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct Args {
  pub width: usize,
  pub height: usize,
  pub zoom: f64,
  pub zpx: f64,
  pub zpy: f64,
  pub iter: u32,
  pub filename: String,
  pub color_map: ColorMap1d,
  pub exponent: f64,
  pub sample_count: u32,
  pub sampler: Sampler,
  #[serde(default)]
  pub post_processing: Vec<PostProcessing>,
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "type")]
pub enum Sampler {
  Uniform,
  /// Kernel density estimation with a uniform kernel.
  ///
  Kde {
    p_min: f64,
    h: f64,
    population: u32,
  },
}

impl Sampler {
  pub fn create_executor(
    self,
    iter: u32,
    viewport: &Viewport,
  ) -> Box<dyn Sampler_<Output = Complex64> + Sync + Send> {
    match self {
      Self::Uniform => Box::new(Uniform::<Complex64>::new()),
      Self::Kde {
        p_min,
        h,
        population,
      } => {
        println!("initializing sampler population");

        let samples = samples(population, p_min, iter, viewport);

        println!("\ninitializing sampler population done");

        let uniform_kde = move |c: &Complex64| {
          let re = (random::<f64>() - 0.5) * h;
          let im = (random::<f64>() - 0.5) * h;

          Complex64::new(c.re + re, c.im + im)
        };

        Box::new(KDE::new(samples, uniform_kde))
      }
    }
  }
}

pub fn buddhabrot(args: Args) {
  let num_pixel = args.width * args.height;

  let buffer = vec_no_clone![AtomicU64::new(0); num_pixel];

  let (w, h) = (args.width as f64, args.height as f64);

  let aspect_ratio = w / h;

  let vp_height = 1. / args.zoom;
  let vp_width = vp_height * aspect_ratio;

  let vp_height_half = vp_height * 0.5;
  let vp_width_half = vp_width * 0.5;

  let x_min = args.zpx - vp_width_half;
  let x_max = vp_width - vp_width_half + args.zpx;

  let y_min = args.zpy - vp_height_half;
  let y_max = vp_height - vp_height_half + args.zpy;

  let viewport = Viewport::new(x_min, y_min, x_max, y_max);

  let delta_x = vp_width / w;
  let delta_y = vp_height / h;

  let sampler = args.sampler.create_executor(args.iter, &viewport);

  println!("starting buddhabrot generation");

  let processed_samples = AtomicU32::new(0);
  (0..args.sample_count).into_par_iter().for_each(|_| {
    let c = sampler.sample();

    let (j, passed_viewport) =
      iter_mandel_check_vp(c, args.iter, &viewport);

    if j != args.iter && passed_viewport {
      let mut z = c;

      for _ in 0..=j {
        let idx = grid_pos(z.re, z.im, delta_x, delta_y, &viewport);

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
  )
  .unwrap();

  println!("successfully written: {}", args.filename);
}

fn samples(
  sample_count: u32,
  p_min: f64,
  iter: u32,
  viewport: &Viewport,
) -> Vec<(Complex64, f64)> {
  let sampler = Uniform::<Complex64>::new();

  let processed_samples = AtomicU32::new(0);

  (0..sample_count)
    .into_par_iter()
    .fold(
      || Vec::new(),
      |mut acc, _| {
        let c = sampler.sample();

        let (j, passed_viewport) =
          iter_mandel_check_vp(c, iter, viewport);

        if j != iter && passed_viewport {
          let p = j as f64 / iter as f64;

          if p > p_min {
            acc.push((c, p));
          }
        }

        let ps = processed_samples.fetch_add(1, Ordering::SeqCst);
        print_progress(ps, sample_count, 2500);

        acc
      },
    )
    .reduce(
      || Vec::new(),
      |mut acc, mut v| {
        acc.append(&mut v);
        acc
      },
    )
}

fn iter_mandel_check_vp(
  c: Complex64,
  iter: u32,
  viewport: &Viewport,
) -> (u32, bool) {
  let mut z = c;
  let mut z_sqr = z.norm_sqr();

  let mut j = 0;
  let mut passed_viewport = false;

  while j < iter && z_sqr <= 4.0 {
    z = z.powi(2) + c;
    z_sqr = z.norm_sqr();

    if viewport.contains_point(z.re, z.im) {
      passed_viewport = true;
    }

    j = j + 1;
  }

  (j, passed_viewport)
}
