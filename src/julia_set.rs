use image::{Rgba, RgbaImage};

use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefMutIterator,
  ParallelIterator,
};

use algorithmic_art::ColorMap1D;

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
    [0, 0, 0, 255],
    [255, 200, 0, 255],
    [255, 255, 255, 255],
    [0, 0, 255, 255],
    [0, 0, 128, 255],
  ]);

  let c = num_complex::Complex::new(-0.4, 0.6);

  //let zoom = 0.0001;
  //let zpx = 0.14;
  //let zpy = -0.39;

  let zoom = 0.5;
  let zpx = 0.;
  let zpy = 0.;

  let vp = 1. / zoom;

  let iter = 1000;

  for x in 0..imgx {
    for y in 0..imgy {
      let zx = x as f32 / imgx as f32;
      let zx = zx * vp - vp / 2. + zpx;
      let zx = zx * norm_x;

      let zy = y as f32 / imgy as f32;
      let zy = zy * vp - vp / 2. + zpy;
      let zy = zy * norm_y;

      let mut z = num_complex::Complex::new(zx, zy);

      let mut i = 0;
      while i < iter && z.norm() <= 2.0 {
        z = z * z + c;
        i += 1;
      }

      let pixel = imgbuf.get_pixel_mut(x, y);

      let color = color_map.value(i as f32 / iter as f32);

      *pixel = Rgba(color);
    }
  }

  imgbuf.save("fractal.png").unwrap();
}
