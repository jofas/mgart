use image::{Rgba, RgbaImage};

use algorithmic_art::{ColorMap1D, RgbaColor};

fn main() {
  let imgx = 1000;
  let imgy = 100;

  let color_map = ColorMap1D::new(vec![
    RgbaColor::new_hex(0xFFFFFFFF),
    RgbaColor::new_hex(0xFF0000FF),
    RgbaColor::new_hex(0x00FF00FF),
    RgbaColor::new_hex(0x0000FFFF),
    RgbaColor::new_hex(0x000000FF),
  ]);

  let mut imgbuf = RgbaImage::new(imgx, imgy);

  for x in 0..imgx {
    for y in 0..imgy {
      let pixel = imgbuf.get_pixel_mut(x, y);
      let color = color_map.value(x as f32 / imgx as f32);
      *pixel = Rgba(color.as_vec());
    }
  }

  imgbuf.save("color_map.png").unwrap();
}
