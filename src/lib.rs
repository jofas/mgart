use image::{Rgba, RgbaImage};

pub mod args;
pub mod util;

use args::JuliaSetArgs;

pub fn julia_set(args: JuliaSetArgs) {
  let mut imgbuf = RgbaImage::new(args.width, args.height);

  let (norm_x, norm_y) = if args.width > args.height {
    (args.width as f32 / args.height as f32, 1.)
  } else {
    (1., args.height as f32 / args.width as f32)
  };

  let vp = 1. / args.zoom;

  for x in 0..args.width {
    for y in 0..args.height {
      let zx = x as f32 / args.width as f32;
      let zx = zx * vp - vp / 2. + args.zpx;
      let zx = zx * norm_x;

      let zy = y as f32 / args.height as f32;
      let zy = zy * vp - vp / 2. + args.zpy;
      let zy = zy * norm_y;

      let mut z = num_complex::Complex::new(zx, zy);

      // if no complex number is given, compute the mandelbrot set
      let c = if let Some(c) = &args.c { c.into() } else { z };

      let mut color = (-z.norm()).exp();

      let mut i = 0;
      while i < args.iter && z.norm() <= 2.0 {
        z = z * z + c;
        color += (-z.norm()).exp();
        i += 1;
      }

      let color = args.color_map.value(color / args.iter as f32);

      let pixel = imgbuf.get_pixel_mut(x, y);
      *pixel = Rgba(color.as_vec());
    }
  }

  imgbuf.save(args.filename).unwrap();
}
