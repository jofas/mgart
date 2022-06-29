use rayon::iter::{
  IndexedParallelIterator, IntoParallelIterator,
  IntoParallelRefMutIterator, ParallelIterator,
};
use rayon::slice::ParallelSliceMut;
use rayon::{current_num_threads, current_thread_index};

use num_complex::Complex64;

use rand::random;

use map_macro::vec_no_clone;

use std::cell::RefCell;
use std::f64::consts::PI;
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, RwLock};

pub mod args;
pub mod util;

use args::{BuddhabrotArgs, ColorMap1dArgs, JuliaSetArgs};

use util::colors::{LCH, RGB};
use util::{Gradient, Sampling, Viewport};

fn attractor(
  z0: Complex64,
  c: Complex64,
  p: usize,
) -> Option<(Complex64, Complex64)> {
  let mut zz = z0;

  for _ in 0..64 {
    let mut z = zz;
    let mut dz = Complex64::new(1., 0.);

    for _ in 0..p {
      dz = 2. * z * dz;
      z = z.powi(2) + c;
    }

    let zz_new = zz - (z - zz) / (dz - 1.);

    if (zz_new - zz).norm_sqr() <= 1e-20 {
      return Some((z, dz));
    }

    zz = zz_new;
  }

  None
}

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

fn samples(
  sample_count: u32,
  p_min: f64,
  iter: u32,
  viewport: &Viewport,
) -> Vec<(Complex64, f64)> {
  let processed_samples = AtomicU32::new(0);
  (0..sample_count)
    .into_par_iter()
    .fold(
      || Vec::new(),
      |mut acc, i| {
        let c = util::random_complex();

        let (j, passed_viewport) =
          iter_mandel_check_vp(c, iter, viewport);

        if j != iter && passed_viewport {
          let p = j as f64 / iter as f64;

          if p > p_min {
            acc.push((c, p));
          }
        }

        let ps = processed_samples.fetch_add(1, Ordering::SeqCst);
        util::print_progress(ps, sample_count, 2500);

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

pub fn buddhabrot(args: BuddhabrotArgs) {
  let num_pixel = args.width * args.height;

  let buffers = vec_no_clone![
    Arc::new(Mutex::new(vec![0.; num_pixel]));
    current_num_threads() * args.buffers_per_thread
  ];

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

  println!("initializing sampler");

  let samples = samples(
    args.sampler.population,
    args.sampler.p_min,
    args.iter,
    &viewport,
  );

  let uniform_kde = |c: &Complex64| {
    let re = (random::<f64>() - 0.5) * args.sampler.h;
    let im = (random::<f64>() - 0.5) * args.sampler.h;

    Complex64::new(c.re + re, c.im + im)
  };

  let sampler = util::KDE::new(samples, uniform_kde);

  println!("\ninitializing sampler done");

  println!("starting buddhabrot generation");

  let processed_samples = AtomicU32::new(0);
  (0..args.sample_count).into_par_iter().for_each(|_| {
    let c = sampler.sample();

    let (j, passed_viewport) =
      iter_mandel_check_vp(c, args.iter, &viewport);

    if j != args.iter && passed_viewport {
      let mut z = c;

      for _ in 0..=j {
        let idx =
          util::grid_pos(z.re, z.im, delta_x, delta_y, &viewport);

        if let Some((x, y)) = idx {
          let b =
            current_thread_index().unwrap() * args.buffers_per_thread;

          let b = b
            + (random::<f64>() * args.buffers_per_thread as f64)
              as usize;

          buffers[b].lock().unwrap()[y * args.width + x] += 1.;
        }

        z = z.powi(2) + c;
      }
    }

    let ps = processed_samples.fetch_add(1, Ordering::SeqCst);
    util::print_progress(ps, args.sample_count, 2500);
  });

  println!("\nbuddhabrot generation finished");

  println!("starting post processing");

  let buffers: Vec<Vec<f64>> = buffers
    .into_iter()
    .map(|b| Arc::try_unwrap(b).unwrap().into_inner().unwrap())
    .collect();

  let (mut min, mut max) = (f64::MAX, 0_f64);

  for b in &buffers {
    let (bmin, bmax) = util::min_max(&buffers[0]);

    min = min.min(bmin);
    max = max.max(bmax);
  }

  let n = buffers.len() as f64;

  let mut avg: Vec<f64> = vec![0.; num_pixel];

  for i in 0..num_pixel {
    let mut sum = 0.;
    let mut squared_sum = 0.;

    for b in &buffers {
      let norm = (b[i] - min) / (max - min);

      sum += norm;
      squared_sum += norm.powi(2);
    }

    avg[i] = sum / n;
  }

  println!("applying gamma correction and color gradient");

  for p in &mut avg {
    *p = args.color_map.gradient().apply_to(*p).powf(args.gamma);
  }

  println!("gamma correction and color gradient applied");

  let buffer = if let Some(smoothing) = args.smoothing {
    println!("starting smoothing process");

    let res = smoothing.smooth(&avg, args.width, args.height);

    println!("\nsmoothing process done");

    res
  } else {
    avg
  };

  println!("post processing done");

  println!("generating final color values");

  let color_map = args.color_map.with_gradient(Gradient::default());

  let mut pixels = vec![0_u8; num_pixel * 3];

  pixels
    .chunks_exact_mut(3)
    .enumerate()
    .for_each(|(i, pixel)| {
      //let count = buffer[i] as f64;
      //let count_norm = (count - min) / (max - min);

      let rgb = color_map.color(buffer[i]);

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

  println!("\nsuccessfully written: {}", args.filename);
}

pub fn julia_set_interior_distance(args: JuliaSetArgs) {
  let num_pixel = args.width * args.height;

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
      let x = (i % args.width) as f64 / w;
      let x = x * vp_width - vp_width_half + args.zpx;

      let y = (i / args.width) as f64 / h;
      let y = y * vp_height - vp_height_half + args.zpy;

      let mut z = Complex64::new(x, y);

      let c = if let Some(c) = &args.c { c.into() } else { z };

      let mut id = None;
      let mut m = f64::MAX;

      for p in 1..=args.iter as usize {
        let z_sqr = z.norm_sqr();

        if z_sqr < m {
          m = z_sqr;

          if let Some((z0, dz0)) = attractor(z, c, p) {
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

        j = j + 1;
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

      let rgb = args.color_map.color(x / args.width as f64);

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
    args.width as u32,
    args.height as u32,
    image::ColorType::Rgb8,
  )
  .unwrap();

  println!("\nsuccessfully written: {}", args.filename);
}
