use image::{Rgba, RgbaImage};

use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefMutIterator,
  ParallelIterator,
};

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

  let c = num_complex::Complex::new(-0.4, 0.6);

  let zoom = 0.01;

  let zpx = 0.14;
  let zpy = -0.4;

  for x in 0..imgx {
    for y in 0..imgy {
      let zx = (x as f32 / imgx as f32 * 2. * zoom - zoom + zpx) * norm_x;
      let zy = (y as f32 / imgy as f32 * 2. * zoom - zoom + zpy) * norm_y;

      let mut z = num_complex::Complex::new(zx, zy);

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
