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

  let mut buf = vec![[0_u8; 3]; num_pixel];

  let pixel_created = Arc::new(AtomicUsize::new(0));

  let (w, h) = (args.width as f64, args.height as f64);

  let aspect_ratio = w / h;

  let vp_height = 1. / args.zoom;
  let vp_width = vp_height * aspect_ratio;

  let vp_height_half = vp_height * 0.5;
  let vp_width_half = vp_width * 0.5;

  let (cx, cy) = if let Some(c) = &args.c {
    (Some(c.re()), Some(c.im()))
  } else {
    (None, None)
  };

  buf.par_iter_mut().enumerate().for_each(|(i, pixel)| {
    let x = (i % args.width) as f64 / w;
    let x = x * vp_width - vp_width_half + args.zpx;

    let y = (i / args.width) as f64 / h;
    let y = y * vp_height - vp_height_half + args.zpy;

    let cx = cx.unwrap_or(x);
    let cy = cy.unwrap_or(y);

    let mut zx = x;
    let mut zy = y;

    let mut zx_sqr = zx.powi(2);
    let mut zy_sqr = zy.powi(2);

    let mut j = 0;
    while j < args.iter && zx_sqr + zy_sqr <= 4.0 {
      zy = (zx * zy) * 2.0 + cy;
      zx = zx_sqr - zy_sqr + cx;

      zx_sqr = zx.powi(2);
      zy_sqr = zy.powi(2);

      j += 1;
    }

    let mu = (zx_sqr + zy_sqr).sqrt().ln().ln() / 2.0_f64.ln();

    let j = (j + 1) as f64 - mu;

    *pixel = args.color_map.color(j / args.iter as f64).as_vec();

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

  let buf: Vec<u8> = buf.into_iter().flatten().collect();

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

  let mut buf = vec![[0_u8; 3]; num_pixel];

  let pixel_created = Arc::new(AtomicUsize::new(0));

  buf.par_iter_mut().enumerate().for_each(|(i, pixel)| {
    let x = i % args.width;

    *pixel =
      args.color_map.color(x as f64 / args.width as f64).as_vec();

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
    image::ColorType::Rgb8,
  )
  .unwrap();

  println!("\nsuccessfully written: {}", args.filename);
}
