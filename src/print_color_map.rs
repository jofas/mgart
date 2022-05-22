use image::{Rgba, RgbaImage};

use algorithmic_art::ColorMap1D;

fn main() {
  let imgx = 1000;
  let imgy = 100;

  let color_map = ColorMap1D::new(vec![
    [255, 255, 255, 255],
    [255, 0, 0, 255],
    [0, 255, 0, 255],
    [0, 0, 255, 255],
    [0, 0, 0, 255],
  ]);

  let mut imgbuf = RgbaImage::new(imgx, imgy);

  for x in 0..imgx {
    for y in 0..imgy {
      let pixel = imgbuf.get_pixel_mut(x, y);
      let color = color_map.value(x as f32 / imgx as f32);
      *pixel = Rgba(color);
    }
  }

  imgbuf.save("color_map.png").unwrap();
}
