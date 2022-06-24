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
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, RwLock};

pub mod args;
pub mod util;

use args::{BuddhabrotArgs, ColorMap1dArgs, JuliaSetArgs};
use util::colors::{LCH, RGB};

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

pub fn buddhabrot(args: BuddhabrotArgs) {
  let num_pixel = args.width * args.height;

  let num_buffers = 8;

  let buffers = vec_no_clone![
    Arc::new(Mutex::new(vec![0.; num_pixel]));
    current_num_threads() * num_buffers
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

  let delta_x = vp_width / w;
  let delta_y = vp_height / h;

  let sampler = RwLock::new(util::GaussianKDE::new(3_000_000, 0.1));

  let processed_samples = AtomicU64::new(0);
  let iter_count = AtomicU64::new(0);
  let hits_count = AtomicU64::new(0);

  (0..args.sample_count).into_par_iter().for_each(|_| {
    let (c, sample_id) = sampler.read().unwrap().sample();
    //let c = util::random_complex();

    let mut z = c;
    let mut z_sqr = z.norm_sqr();

    let mut j = 0;
    let mut passed_viewport = false;
    while j < args.iter && z_sqr <= 4.0 {
      z = z.powi(2) + c;
      z_sqr = z.norm_sqr();

      if util::in_viewport(z.re, z.im, x_min, x_max, y_min, y_max) {
        passed_viewport = true;
      }

      j = j + 1;
    }

    let p = if j != args.iter && passed_viewport {
      iter_count.fetch_add(j as u64, Ordering::SeqCst);
      hits_count.fetch_add(1, Ordering::SeqCst);

      j as f64 / args.iter as f64
    } else {
      0.
    };

    sampler.write().unwrap().update(c, p, sample_id);

    if p > 0. {
      let mut z = c;

      for i in 0..=j {
        let idx = util::grid_pos(
          z.re, z.im, x_min, x_max, y_min, y_max, delta_x, delta_y,
        );

        if let Some((x, y)) = idx {
          let b = current_thread_index().unwrap() * num_buffers;
          let b = b + (random::<f64>() * num_buffers as f64) as usize;

          buffers[b].lock().unwrap()[y * args.width + x] += 1.;
            //p - i as f64 / args.iter as f64;
        }

        z = z.powi(2) + c;
      }
    }

    let ps = processed_samples.fetch_add(1, Ordering::SeqCst);

    if ps % 2500 == 0 || ps == args.sample_count - 1 {
      print!(
        "{}/{} samples iterated ({:.2}%)\r",
        ps + 1,
        args.sample_count,
        (ps as f64 / args.sample_count as f64) * 100.,
      );
    }
  });

  let avg_iter = iter_count.into_inner() / args.iter as u64;
  let avg_iter = avg_iter as f64 / args.sample_count as f64;

  let avg_hits =
    hits_count.into_inner() as f64 / args.sample_count as f64;

  println!("\naverage hits: {}", avg_hits);
  println!("average iterations: {}", avg_iter);
  println!(
    "average probability: {}",
    sampler.into_inner().unwrap().average_probability()
  );

  let buffers: Vec<Vec<f64>> = buffers
    .into_iter()
    .map(|b| Arc::try_unwrap(b).unwrap().into_inner().unwrap())
    .collect();

  let n = buffers.len() as f64;

  let mut avg: Vec<f64> = vec![0.; num_pixel];
  let mut variance: Vec<f64> = vec![0.; num_pixel];

  for i in 0..num_pixel {
    let mut sum = 0.;
    let mut squared_sum = 0.;

    for b in &buffers {
      sum += b[i];
      squared_sum += b[i].powi(2);
    }

    avg[i] = sum / n;
    variance[i] = squared_sum / n - (sum / n).powi(2);
  }

  let mut smoothed: Vec<f64> = vec![0.; num_pixel];

  let mut sat: Vec<f64> = vec![avg[0]; num_pixel];

  for i in 1..num_pixel {
    let x = i % args.width;
    let y = i / args.width;

    // the value to the left
    sat[i] += sat[i - 1];

    // the value above
    if y > 0 {
      sat[i] += sat[(y - 1) * args.width + x];
    }
  }

  println!("starting smoothing process");

  /*
  let n = 11;
  let window_size = 51;

  let processed_pixels = AtomicUsize::new(0);

  smoothed.par_iter_mut().enumerate().for_each(|(i, pixel)| {
    let mut s = 0.;
    let mut cp = 0.;

    let x = i % args.width;
    let y = i / args.width;

    let (wx0, wy0, wx1, wy1) = util::discrete_bounded_square(
      x,
      y,
      window_size,
      args.width - 1,
      args.height - 1,
    );

    let (x0, y0, x1, y1) = util::discrete_bounded_square(
      x,
      y,
      n,
      args.width - 1,
      args.height - 1,
    );

    let bp = sat[y1 * args.width + x1] + sat[y0 * args.width + x0]
      - sat[y1 * args.width + x0] - sat[y0 * args.width + x1];
    let bp = bp / (x1 - x0 + 1) as f64 / (y1 - y0 + 1) as f64;

    for x in wx0..=wx1 {
      for y in wy0..=wy1 {
        let (x0, y0, x1, y1) = util::discrete_bounded_square(
          x,
          y,
          n,
          args.width - 1,
          args.height - 1,
        );

        let bq = sat[y1 * args.width + x1] + sat[y0 * args.width + x0]
          - sat[y1 * args.width + x0] - sat[y0 * args.width + x1];
        let bq = bq / (x1 - x0 + 1) as f64 / (y1 - y0 + 1) as f64;

        let fpq = (-((bq - bp).powi(2) / variance[i])).exp();

        s += avg[i] * fpq;
        cp += fpq;
      }
    }

    *pixel = s / cp;

    let pc = processed_pixels.fetch_add(1, Ordering::SeqCst);

    if pc % 10 == 0 || pc == num_pixel - 1 {
      print!(
        "{}/{} pixels smoothed ({:.2}%)\r",
        pc + 1,
        num_pixel,
        (pc as f64 / num_pixel as f64) * 100.,
      );
    }
  });
  */

  let counter_buf = avg; //smoothed;

  let mut max = 0.;
  let mut min = f64::MAX;

  for counter in &counter_buf {
    let counter = *counter;

    if counter < min {
      min = counter;
    }

    if counter > max {
      max = counter;
    }
  }

  let mut buf = vec![0_u8; num_pixel * 3];

  buf.chunks_exact_mut(3).enumerate().for_each(|(i, pixel)| {
    let count = counter_buf[i] as f64;
    let count_norm = (count - min) / (max - min);

    let rgb = args.color_map.color(count_norm);

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
