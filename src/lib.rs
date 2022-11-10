#![allow(clippy::must_use_candidate)]

use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefMutIterator,
  ParallelIterator,
};
use rayon::slice::ParallelSliceMut;

use serde::Deserialize;

use num_complex::Complex64;

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

pub mod args;
pub mod buddhabrot;
pub mod util;

use crate::args::{ColorMap1dArgs, JuliaSetArgs};

use crate::util::coloring::colors::RGB;
use crate::util::finite_attractor;

use crate::buddhabrot::{buddhabrot, Args as BuddhabrotArgs};

#[derive(Deserialize)]
#[serde(tag = "algorithm")]
#[serde(rename_all = "snake_case")]
pub enum Algorithm {
  JuliaSet(JuliaSetArgs),
  Buddhabrot(BuddhabrotArgs),
  #[serde(rename = "color_map_1d")]
  ColorMap1d(ColorMap1dArgs),
}

impl Algorithm {
  pub fn execute(self) {
    match self {
      Self::JuliaSet(args) => {
        println!("generating julia set with arguments:\n{}", args);
        julia_set_interior_distance(args);
      }
      Self::Buddhabrot(args) => {
        println!("generating buddhabrot with arguments: \n{}", args);
        buddhabrot(args);
      }
      Self::ColorMap1d(args) => {
        println!("generating 1d color map with arguments:\n{}", args);
        color_map_1d(args);
      }
    }
  }
}

#[derive(Deserialize)]
pub struct Algorithms(Vec<Algorithm>);

impl Algorithms {
  pub fn execute(self) {
    for cmd in self.0 {
      cmd.execute();
    }
  }
}

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

pub fn julia_set_interior_distance(args: JuliaSetArgs) {
  let (width, height) = (args.width as usize, args.height as usize);

  let num_pixel = width * height;

  let mut distance_buf = vec![0.; num_pixel];

  let d_max = Arc::new(Mutex::new(f64::MIN));

  let (w, h) = (args.width as f64, args.height as f64);

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

        let mut d_max = d_max.lock().unwrap();

        if id > *d_max {
          *d_max = id;
        }
      } else {
        *pixel = 0.;
      }
    });

  let mut buf = vec![0_u8; num_pixel * 3];

  let d_max = Arc::try_unwrap(d_max).unwrap().into_inner().unwrap();

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
  )
  .unwrap();

  println!("\nsuccessfully written: {}", args.filename);
}

pub fn julia_set(args: JuliaSetArgs) {
  let (width, height) = (args.width as usize, args.height as usize);

  let num_pixel = width * height;

  let mut buf = vec![0_u8; num_pixel * 3];

  let pixel_created = Arc::new(AtomicUsize::new(0));

  let (w, h) = (f64::from(args.width), f64::from(args.height));

  let aspect_ratio = w / h;

  let vp_height = 1. / args.zoom;
  let vp_width = vp_height * aspect_ratio;

  let vp_height_half = vp_height * 0.5;
  let vp_width_half = vp_width * 0.5;

  buf
    .par_chunks_exact_mut(3)
    .enumerate()
    .for_each(|(i, pixel)| {
      let x = (i % width) as f64 / w;
      let x = x * vp_width - vp_width_half + args.zpx;

      let y = (i / width) as f64 / h;
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

      //let rgb = args.color_map.color(color);

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
    args.width,
    args.height,
    image::ColorType::Rgb8,
  )
  .unwrap();

  println!("\nsuccessfully written: {}", args.filename);
}

pub fn color_map_1d(args: ColorMap1dArgs) {
  let (w, h) = (args.width as usize, args.height as usize);

  let num_pixel = w * h;

  let mut buf = vec![0_u8; num_pixel * 3];

  let pixel_created = Arc::new(AtomicUsize::new(0));

  buf
    .par_chunks_exact_mut(3)
    .enumerate()
    .for_each(|(i, pixel)| {
      let x = (i % w) as f64;

      let rgb = args.color_map.color(x / w as f64);

      pixel[0] = rgb.r();
      pixel[1] = rgb.g();
      pixel[2] = rgb.b();

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
    args.width,
    args.height,
    image::ColorType::Rgb8,
  )
  .unwrap();

  println!("\nsuccessfully written: {}", args.filename);
}
