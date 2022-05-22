use image::codecs::gif::{GifEncoder, Repeat};
use image::{Frame, Rgba, RgbaImage};

use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefMutIterator,
  ParallelIterator,
};

use std::f32::consts::PI;
use std::fs::File;

fn create_frame(frame: &mut Frame, a: f32) {
  let imgbuf = frame.buffer_mut();

  let (imgx, imgy) = imgbuf.dimensions();

  let scalex = 3.0 / imgx as f32;
  let scaley = 3.0 / imgy as f32;

  for x in 0..imgx {
    for y in 0..imgy {
      let cx = x as f32 * scalex - 1.5;
      let cy = y as f32 * scaley - 1.5;

      //println!("x: {}, y: {}, cx: {}, cy: {}", x, y, cx, cy);

      //let c = num_complex::Complex::new(-0.4, 0.6);
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
}

fn main() {
  let imgx = 1080;
  let imgy = 1080;

  let num_frames = 1_000;

  let mut frames =
    vec![Frame::new(RgbaImage::new(imgx, imgy)); num_frames];

  frames.par_iter_mut().enumerate().for_each(|(i, f)| {
    println!("creating frame: {}", i);
    let a = i as f32 * 2. * PI / num_frames as f32;
    create_frame(f, a);
  });

  let mut gif = GifEncoder::new_with_speed(
    File::create("fractal.gif").unwrap(),
    1,
  );

  gif.set_repeat(Repeat::Infinite).unwrap();

  println!("writing frames to file");

  gif.encode_frames(frames).unwrap();

  println!("finished creating gif");
}
