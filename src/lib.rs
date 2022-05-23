use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefMutIterator,
  ParallelIterator,
};

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

pub mod args;
pub mod util;

use args::JuliaSetArgs;

pub fn julia_set(args: JuliaSetArgs) {
  let num_pixel = args.width * args.height;

  let mut buf = vec![[0_u8; 4]; num_pixel];

  let pixel_created = Arc::new(AtomicUsize::new(0));

  let (norm_x, norm_y) = if args.width > args.height {
    (args.width as f64 / args.height as f64, 1.)
  } else {
    (1., args.height as f64 / args.width as f64)
  };

  let vp = 1. / args.zoom;

  buf.par_iter_mut().enumerate().for_each(|(i, pixel)| {
    let x = i % args.width;
    let y = i / args.width;

    let zx = x as f64 / args.width as f64;
    let zx = zx * vp - vp / 2. + args.zpx;
    let zx = zx * norm_x;

    let zy = y as f64 / args.height as f64;
    let zy = zy * vp - vp / 2. + args.zpy;
    let zy = zy * norm_y;

    let mut z = num_complex::Complex::new(zx, zy);

    // if no complex number is given, compute the mandelbrot set
    let c = if let Some(c) = &args.c { c.into() } else { z };

    let mut color = (-z.norm()).exp();

    let mut j = 0;
    while j < args.iter && z.norm() <= 2.0 {
      z = z * z + c;
      color += (-z.norm()).exp();
      j += 1;
    }

    *pixel = args.color_map.value(color / args.iter as f64).as_vec();

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
    args.filename,
    &buf,
    args.width as u32,
    args.height as u32,
    image::ColorType::Rgba8,
  )
  .unwrap()
}
