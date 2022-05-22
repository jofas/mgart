use image::{Rgba, RgbaImage};

use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefMutIterator,
  ParallelIterator,
};

use algorithmic_art::{ColorMap1D, RgbaColor};

fn main() {
  //let imgx = 7680;
  //let imgy = 4320;

  let imgx = 1920;
  let imgy = 1080;

  let mut imgbuf = RgbaImage::new(imgx, imgy);

  let (norm_x, norm_y) = if imgx > imgy {
    (imgx as f32 / imgy as f32, 1.)
  } else {
    (1., imgy as f32 / imgx as f32)
  };

  let color_map = ColorMap1D::new(vec![
    RgbaColor::new_hex(0x000000FF),
    RgbaColor::new_hex(0x14213DFF),
    RgbaColor::new_hex(0xFCA311FF),
    RgbaColor::new_hex(0xE5E5E5FF),
    RgbaColor::new_hex(0xFFFFFFFF),
  ]);

  let c = num_complex::Complex::new(-0.4, 0.6);

  let zoom = 0.5;
  let zpx = 0.;
  let zpy = 0.;

  let zoom = 150.;
  let zpx = 0.138;
  let zpy = 0.142;

  let vp = 1. / zoom;

  let iter = 250;

  for x in 0..imgx {
    for y in 0..imgy {
      let zx = x as f32 / imgx as f32;
      let zx = zx * vp - vp / 2. + zpx;
      let zx = zx * norm_x;

      let zy = y as f32 / imgy as f32;
      let zy = zy * vp - vp / 2. + zpy;
      let zy = zy * norm_y;

      let mut z = num_complex::Complex::new(zx, zy);
      let mut color = (-z.norm()).exp();

      let mut i = 0;
      while i < iter && z.norm() <= 2.0 {
        z = z * z + c;
        color += (-z.norm()).exp();
        i += 1;
      }

      let color = color_map.value(color / iter as f32);

      let pixel = imgbuf.get_pixel_mut(x, y);
      *pixel = Rgba(color.as_vec());
    }
  }

  imgbuf.save("fractal.png").unwrap();
}
