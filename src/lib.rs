use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefMutIterator,
  ParallelIterator,
};

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

pub mod args;
pub mod util;

use args::{ColorMap1dArgs, JuliaSetArgs};

pub fn julia_set(args: JuliaSetArgs) {
  let num_pixel = args.width * args.height;

  let mut buf = vec![[0_u8; 4]; num_pixel];

  let pixel_created = Arc::new(AtomicUsize::new(0));

  let (w, h) = (args.width as f64, args.height as f64);

  let aspect_ratio = w / h;

  let vp = 1. / args.zoom;

  buf.par_iter_mut().enumerate().for_each(|(i, pixel)| {
    // instead of doing this as float, do it as int and only
    // afterwards scale down
    //
    let zx = (i % args.width) as f64 / w;
    let zx = zx * vp - vp / 2. + args.zpx;
    //let zx = zx * aspect_ratio;

    let zy = (i / args.width) as f64 / h;
    let zy = zy * vp - vp / 2. + args.zpy;

    let mut z = num_complex::Complex::new(zx, zy);

    // if no complex number is given, compute the mandelbrot set
    let c = if let Some(c) = &args.c { c.into() } else { z };

    let mut j = 0;
    while j < args.iter && z.norm_sqr() <= 4.0 {
      z = z * z + c;
      j += 1;
    }

    let nu =
      ((z.norm_sqr().ln() / 2.) / 2.0_f64.ln()).ln() / 2.0_f64.ln();

    let j = j as f64 + 1. - nu;

    *pixel = args
      .color_map
      .color(j / args.iter as f64, &args.color_method)
      .as_vec();

    let pc = pixel_created.fetch_add(1, Ordering::SeqCst);

    print!(
      "{}/{} pixels created ({:.2}%)\r",
      pc,
      num_pixel,
      (pc as f32 / num_pixel as f32) * 100.,
    );
  });

  let buf: Vec<u8> = buf.into_iter().flatten().collect();

  image::save_buffer(
    &args.filename,
    &buf,
    args.width as u32,
    args.height as u32,
    image::ColorType::Rgba8,
  )
  .unwrap();

  println!("\nsuccessfully written: {}", args.filename);
}

pub fn color_map_1d(args: ColorMap1dArgs) {
  let num_pixel = args.width * args.height;

  let mut buf = vec![[0_u8; 4]; num_pixel];

  let pixel_created = Arc::new(AtomicUsize::new(0));

  buf.par_iter_mut().enumerate().for_each(|(i, pixel)| {
    let x = i % args.width;

    *pixel = args
      .color_map
      .color(x as f64 / args.width as f64, &args.color_method)
      .as_vec();

    let pc = pixel_created.fetch_add(1, Ordering::SeqCst);

    print!(
      "{}/{} pixels created ({:.2}%)\r",
      pc,
      num_pixel,
      (pc as f32 / num_pixel as f32) * 100.,
    );
  });

  let buf: Vec<u8> = buf.into_iter().flatten().collect();

  image::save_buffer(
    &args.filename,
    &buf,
    args.width as u32,
    args.height as u32,
    image::ColorType::Rgba8,
  )
  .unwrap();

  println!("\nsuccessfully written: {}", args.filename);
}
