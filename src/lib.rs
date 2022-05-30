use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use rayon::slice::ParallelSliceMut;

use num_complex::Complex64;

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

pub mod args;
pub mod util;

use args::{ColorMap1dArgs, JuliaSetArgs};
use util::colors::{LCH, RGB};

fn attractor(
  zx: f64,
  zy: f64,
  cx: f64,
  cy: f64,
  p: usize,
) -> Option<(f64, f64, f64, f64)> {
  unimplemented!()
}

pub fn julia_set(args: JuliaSetArgs) {
  let num_pixel = args.width * args.height;

  let mut buf = vec![0_u8; num_pixel * 3];

  let pixel_created = Arc::new(AtomicUsize::new(0));

  let (w, h) = (args.width as f64, args.height as f64);

  let aspect_ratio = w / h;

  let vp_height = 1. / args.zoom;
  let vp_width = vp_height * aspect_ratio;

  let vp_height_half = vp_height * 0.5;
  let vp_width_half = vp_width * 0.5;

  buf
    .par_chunks_exact_mut(3)
    .enumerate()
    .for_each(|(i, pixel)| {
      let x = (i % args.width) as f64 / w;
      let x = x * vp_width - vp_width_half + args.zpx;

      let y = (i / args.width) as f64 / h;
      let y = y * vp_height - vp_height_half + args.zpy;

      let mut z = Complex64::new(x, y);

      let c = if let Some(c) = &args.c { c.into() } else { z };

      let mut z_sqr = z.norm_sqr();

      // outer distance
      //let mut dzx = 0.;
      //let mut dzy = 0.;

      let mut j = 0;
      while j < args.iter && z_sqr <= 4.0 {
        //dzx = 2. * zx * dzx + 1.;
        //dzy = 2. * zy * dzy + 1.;

        z = z.powi(2) + c;
        z_sqr = z.norm_sqr();
        j += 1;
      }

      let color = if j == args.iter {
        1.
        //j as f64 / args.iter as f64
      } else {
        let mu = z_sqr.sqrt().log2().log2();
        ((j + 1) as f64 - mu) / args.iter as f64

        /*
        let z_mag = (zx_sqr + zy_sqr).sqrt();
        let dz_mag = (dzx.powi(2) + dzy.powi(2)).sqrt();
        //let distance = z_mag.powi(2).ln() * z_mag / dz_mag;
        //let distance = 0. - 5. * distance.ln() / args.zoom.ln();

        let distance = 2. * z_mag * z_mag.ln() / dz_mag;

        //println!("distance: {}", distance);
        distance
        */
      };

      let rgb = args.color_map.color(color).as_vec();

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

      let pc = pixel_created.fetch_add(1, Ordering::SeqCst);

      if pc % 2500 == 0 || pc == num_pixel {
        print!(
          "{}/{} pixels created ({:.2}%)\r",
          pc,
          num_pixel,
          (pc as f32 / num_pixel as f32) * 100.,
        );
      }
    });

  image::save_buffer(
    &args.filename,
    &buf,
    args.width as u32,
    args.height as u32,
    image::ColorType::Rgb8,
  )
  .unwrap();

  println!("\nsuccessfully written: {}", args.filename);
}

pub fn color_map_1d(args: ColorMap1dArgs) {
  let num_pixel = args.width * args.height;

  let mut buf = vec![0_u8; num_pixel * 3];

  let pixel_created = Arc::new(AtomicUsize::new(0));

  buf
    .par_chunks_exact_mut(3)
    .enumerate()
    .for_each(|(i, pixel)| {
      let x = (i % args.width) as f64;

      let rgb = args.color_map.color(x / args.width as f64).as_vec();

      pixel[0] = rgb[0];
      pixel[1] = rgb[1];
      pixel[2] = rgb[2];

      let pc = pixel_created.fetch_add(1, Ordering::SeqCst);

      print!(
        "{}/{} pixels created ({:.2}%)\r",
        pc,
        num_pixel,
        (pc as f32 / num_pixel as f32) * 100.,
      );
    });

  image::save_buffer(
    &args.filename,
    &buf,
    args.width as u32,
    args.height as u32,
    image::ColorType::Rgb8,
  )
  .unwrap();

  println!("\nsuccessfully written: {}", args.filename);
}
