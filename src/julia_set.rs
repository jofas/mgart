use serde::{Deserialize, Serialize};

use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use rayon::slice::ParallelSliceMut;

use display_json::DisplayAsJsonPretty;

use num_complex::Complex64;

use anyhow::Result;

use log::info;

use crate::util::coloring::colors::RGB;
use crate::util::coloring::ColorMap1d;
use crate::util::{ComplexNumber, ProgressPrinter};

/*
fn interior_distance(z0: Complex64, c: Complex64, p: usize) -> f64 {
  let mut z = z0;
  let mut dz = Complex64::new(1., 0.);
  let mut dzdz = Complex64::new(0., 0.);
  let mut dc = Complex64::new(0., 0.);
  let mut dcdz = Complex64::new(0., 0.);

  for _ in 0..p {
    dcdz = 2. * (z * dcdz + dz * dc);
    dc = 2. * z * dc + 1.;
    dzdz = 2. * (dz * dz + z * dzdz);
    dz = 2. * z * dz;
    z = z.powi(2) + c;
  }

  (1. - dz.norm_sqr()) / (dcdz + dzdz * dc / (1. - dz)).norm()
}
*/

/*
/// Creates a rendering of a julia set as a `PNG` image.
///
/// # Errors
///
/// Returns an error, if the generated `PNG` image could not be saved
/// to disk.
///
pub fn julia_set_interior_distance(
  args: &JuliaSetArgs,
) -> Result<()> {
  let (width, height) = (args.width as usize, args.height as usize);

  let num_pixel = width * height;

  let mut distance_buf = vec![0.; num_pixel];

  let d_max = Arc::new(Mutex::new(f64::MIN));

  let (w, h) = (f64::from(args.width), f64::from(args.height));

  let aspect_ratio = w / h;

  let vp_height = 1. / args.zoom;
  let vp_width = vp_height * aspect_ratio;

  let vp_height_half = vp_height * 0.5;
  let vp_width_half = vp_width * 0.5;

  distance_buf
    .par_iter_mut()
    .enumerate()
    .for_each(|(i, pixel)| {
      let x = (i % width) as f64 / w;
      let x = x * vp_width - vp_width_half + args.zpx;

      let y = (i / width) as f64 / h;
      let y = y * vp_height - vp_height_half + args.zpy;

      let mut z = Complex64::new(x, y);

      let c = if let Some(c) = &args.c { c.into() } else { z };

      let mut id = None;
      let mut m = f64::MAX;

      for p in 1..=args.iter as usize {
        let z_sqr = z.norm_sqr();

        if z_sqr < m {
          m = z_sqr;

          if let Some((_z0, dz0)) = finite_attractor(z, c, p) {
            if dz0.norm_sqr() <= 1. {
              id = Some(dz0.norm_sqr());
              // dz0 = inner coordinate
              // z0 = finite attractor
              //id = Some(interior_distance(z0, c, p));
              //id = Some(p as f64);
              break;
            }
          }
        }

        if z_sqr >= 4.0 {
          break;
        }

        z = z.powi(2) + c;
      }

      if let Some(id) = id {
        let id = id.abs();

        *pixel = id;

        if let Ok(mut d_max) = d_max.lock() {
          if id > *d_max {
            *d_max = id;
          }
        }
      } else {
        *pixel = 0.;
      }
    });

  let mut buf = vec![0_u8; num_pixel * 3];

  let d_max = Arc::try_unwrap(d_max)
    .map_err(|_| anyhow!("impossible"))?
    .into_inner()?;

  buf
    .par_chunks_exact_mut(3)
    .enumerate()
    .for_each(|(i, pixel)| {
      let d = distance_buf[i];
      let d_norm = (d_max - d) / d_max;

      let rgb = args.color_map.color(d_norm);

      pixel[0] = rgb.r();
      pixel[1] = rgb.g();
      pixel[2] = rgb.b();
    });

  image::save_buffer(
    &args.filename,
    &buf,
    args.width as u32,
    args.height as u32,
    image::ColorType::Rgb8,
  )?;

  info!("\nsuccessfully written: {}", args.filename);

  Ok(())
}
*/

#[derive(Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct JuliaSet {
  pub width: u32,
  pub height: u32,
  pub zoom: f64,
  pub zpx: f64,
  pub zpy: f64,
  pub iter: u64,
  pub filename: String,
  pub color_map: ColorMap1d,
  pub c: Option<ComplexNumber>,
}

impl JuliaSet {
  /// Creates a rendering of a julia set as a `PNG` image.
  ///
  /// # Errors
  ///
  /// Returns an error, if the generated `PNG` image could not be saved
  /// to disk.
  ///
  pub fn create(&self) -> Result<()> {
    let (width, height) = (self.width as usize, self.height as usize);

    let num_pixel = width * height;

    let mut buf = vec![0_u8; num_pixel * 3];

    let pp = ProgressPrinter::new(num_pixel as u64, 2500);

    let (w, h) = (f64::from(self.width), f64::from(self.height));

    let aspect_ratio = w / h;

    let vp_height = 1. / self.zoom;
    let vp_width = vp_height * aspect_ratio;

    let vp_height_half = vp_height * 0.5;
    let vp_width_half = vp_width * 0.5;

    buf
      .par_chunks_exact_mut(3)
      .enumerate()
      .for_each(|(i, pixel)| {
        let x = (i % width) as f64 / w;
        let x = x * vp_width - vp_width_half + self.zpx;

        let y = (i / width) as f64 / h;
        let y = y * vp_height - vp_height_half + self.zpy;

        let mut z = Complex64::new(x, y);

        let c = if let Some(c) = &self.c { c.into() } else { z };

        let mut z_sqr = z.norm_sqr();

        // outer distance
        //let mut dzx = 0.;
        //let mut dzy = 0.;

        let mut j = 0;
        while j < self.iter && z_sqr <= 4.0 {
          //dzx = 2. * zx * dzx + 1.;
          //dzy = 2. * zy * dzy + 1.;

          z = z.powi(2) + c;
          z_sqr = z.norm_sqr();

          j += 1;
        }

        let color = if j == self.iter {
          1.
          //j as f64 / self.iter as f64
        } else {
          let mu = z_sqr.sqrt().log2().log2();
          ((j + 1) as f64 - mu) / self.iter as f64

          /*
          let z_mag = (zx_sqr + zy_sqr).sqrt();
          let dz_mag = (dzx.powi(2) + dzy.powi(2)).sqrt();
          //let distance = z_mag.powi(2).ln() * z_mag / dz_mag;
          //let distance = 0. - 5. * distance.ln() / self.zoom.ln();

          let distance = 2. * z_mag * z_mag.ln() / dz_mag;

          //info!("distance: {}", distance);
          distance
          */
        };

        //let rgb = self.color_map.color(color);

        /*
        let rgb = LCH::new(
          //((color * 100. * std::f64::consts::PI).sin() / 2. + 0.5) * 100.,
          color.powf(3.5).fract() * 100.,
          0.,
          0.,
          //65.,
          //100. - (100. * color),
          //132.,//32. + 100. - (100. * color),
          //(360. * color).powi(2),
        ).rgb();
        */

        let rgb = RGB::new(
          (color * 255.) as u8,
          (color * 255.) as u8,
          (color * 255.) as u8,
        );

        pixel[0] = rgb.r();
        pixel[1] = rgb.g();
        pixel[2] = rgb.b();

        pp.increment();
      });

    image::save_buffer(
      &self.filename,
      &buf,
      self.width,
      self.height,
      image::ColorType::Rgb8,
    )?;

    info!("successfully written: {}", self.filename);

    Ok(())
  }
}
