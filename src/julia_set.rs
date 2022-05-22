use image::{Rgba, RgbaImage};

use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefMutIterator,
  ParallelIterator,
};

fn main() {
  let imgx = 1080;
  let imgy = 1080;

  let mut imgbuf = RgbaImage::new(imgx, imgy);

  let scalex = 3.0 / imgx as f32;
  let scaley = 3.0 / imgy as f32;

  let c = num_complex::Complex::new(-0.4, 0.6);

  for x in 0..imgx {
    for y in 0..imgy {
      let cx = x as f32 * scalex - 1.5;
      let cy = y as f32 * scaley - 1.5;

      let mut z = num_complex::Complex::new(cx, cy);

      let mut i = 0;
      while i < 255 && z.norm() <= 2.0 {
        z = z * z + c;
        i += 1;
      }

      let pixel = imgbuf.get_pixel_mut(x, y);

      if i == 255 {
        *pixel = Rgba([0, 0, 0, 255]);
      } else {
        *pixel = Rgba([0, i as u8, 0, 255]);
      }
    }
  }

  imgbuf.save("fractal.png").unwrap();
}
