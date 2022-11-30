use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use num::cast;
use num::integer::div_rem;

use crate::util::coloring::colors::{Color, RGB};
use crate::util::coloring::ColorMap1d;
use crate::util::frame::Frame;
use crate::util::viewport::Viewport;
use crate::util::{ComplexNumber, ProgressPrinter};

/*
fn interior_distance(z0: Complex64, c: Complex64, p: usize) -> f64 {
  let mut z = z0;
  let mut dz = Complex64::new(1., 0.);
  let mut dzdz = Complex64::new(0., 0.);
  let mut dc = Complex64::new(0., 0.);
  let mut dcdz = Complex64::new(0., 0.);

  for _ in 0..p {
    dcdz = 2. * (z * dcdz + dz * dc);
    dc = 2. * z * dc + 1.;
    dzdz = 2. * (dz * dz + z * dzdz);
    dz = 2. * z * dz;
    z = z.powi(2) + c;
  }

  (1. - dz.norm_sqr()) / (dcdz + dzdz * dc / (1. - dz)).norm()
}
*/

/*
/// Creates a rendering of a julia set as a `PNG` image.
///
/// # Errors
///
/// Returns an error, if the generated `PNG` image could not be saved
/// to disk.
///
pub fn julia_set_interior_distance(
  args: &JuliaSetArgs,
) -> Result<()> {
  let (width, height) = (args.width as usize, args.height as usize);

  let num_pixel = width * height;

  let mut distance_buf = vec![0.; num_pixel];

  let d_max = Arc::new(Mutex::new(f64::MIN));

  let (w, h) = (f64::from(args.width), f64::from(args.height));

  let aspect_ratio = w / h;

  let vp_height = 1. / args.zoom;
  let vp_width = vp_height * aspect_ratio;

  let vp_height_half = vp_height * 0.5;
  let vp_width_half = vp_width * 0.5;

  distance_buf
    .par_iter_mut()
    .enumerate()
    .for_each(|(i, pixel)| {
      let x = (i % width) as f64 / w;
      let x = x * vp_width - vp_width_half + args.zpx;

      let y = (i / width) as f64 / h;
      let y = y * vp_height - vp_height_half + args.zpy;

      let mut z = Complex64::new(x, y);

      let c = if let Some(c) = &args.c { c.into() } else { z };

      let mut id = None;
      let mut m = f64::MAX;

      for p in 1..=args.iter as usize {
        let z_sqr = z.norm_sqr();

        if z_sqr < m {
          m = z_sqr;

          if let Some((_z0, dz0)) = finite_attractor(z, c, p) {
            if dz0.norm_sqr() <= 1. {
              id = Some(dz0.norm_sqr());
              // dz0 = inner coordinate
              // z0 = finite attractor
              //id = Some(interior_distance(z0, c, p));
              //id = Some(p as f64);
              break;
            }
          }
        }

        if z_sqr >= 4.0 {
          break;
        }

        z = z.powi(2) + c;
      }

      if let Some(id) = id {
        let id = id.abs();

        *pixel = id;

        if let Ok(mut d_max) = d_max.lock() {
          if id > *d_max {
            *d_max = id;
          }
        }
      } else {
        *pixel = 0.;
      }
    });

  let mut buf = vec![0_u8; num_pixel * 3];

  let d_max = Arc::try_unwrap(d_max)
    .map_err(|_| anyhow!("impossible"))?
    .into_inner()?;

  buf
    .par_chunks_exact_mut(3)
    .enumerate()
    .for_each(|(i, pixel)| {
      let d = distance_buf[i];
      let d_norm = (d_max - d) / d_max;

      let rgb = args.color_map.color(d_norm);

      pixel[0] = rgb.r();
      pixel[1] = rgb.g();
      pixel[2] = rgb.b();
    });

  image::save_buffer(
    &args.filename,
    &buf,
    args.width as u32,
    args.height as u32,
    image::ColorType::Rgb8,
  )?;

  info!("\nsuccessfully written: {}", args.filename);

  Ok(())
}
*/

#[derive(Serialize, Deserialize, DisplayAsJsonPretty)]
pub struct JuliaSet {
  width: u32,
  height: u32,
  center: ComplexNumber,
  zoom: f64,
  iter: u32,
  rotation: Option<u16>,
  color_map: ColorMap1d,
  c: Option<ComplexNumber>,
}

impl JuliaSet {
  /// Transforms `self` into a [`Creator`] that can be used
  /// to create a rendering of a julia set or the mandelbrot set.
  ///
  #[must_use]
  pub fn creator(self) -> Creator {
    Creator::new(self)
  }
}

pub struct Creator {
  args: JuliaSet,
  viewport: Viewport,
}

impl Creator {
  /// Creates a new instance of [`Creator`].
  ///
  #[must_use]
  pub fn new(args: JuliaSet) -> Self {
    let viewport = Self::viewport(&args);

    Self { args, viewport }
  }

  /// Creates a rendering of a julia set as a `PNG` image.
  ///
  /// # Panics
  ///
  /// Panics, if the configuration is faulty, e.g if there happens an
  /// overflow.
  ///
  #[must_use]
  pub fn create(&self) -> Frame<Color> {
    let mut frame =
      Frame::filled_default(self.args.width, self.args.height);

    let pp = u64::from(self.args.width * self.args.height);
    let pp = ProgressPrinter::new(pp, 2500);

    frame.par_for_each_mut(|(i, pixel)| {
      let (im, re) = div_rem(i, self.args.width.try_into().unwrap());

      let mut z = self.viewport.rotated_point(re, im);

      let c = if let Some(c) = &self.args.c {
        c.into()
      } else {
        z
      };

      let mut z_sqr = z.norm_sqr();

      // outer distance
      //let mut dzx = 0.;
      //let mut dzy = 0.;

      let mut j = 0;
      while j < self.args.iter && z_sqr <= 4.0 {
        //dzx = 2. * zx * dzx + 1.;
        //dzy = 2. * zy * dzy + 1.;

        z = z.powi(2) + c;
        z_sqr = z.norm_sqr();

        j += 1;
      }

      let color = if j == self.args.iter {
        1.
        //j as f64 / self.args.iter as f64
      } else {
        let mu = z_sqr.sqrt().log2().log2();
        (f64::from(j + 1) - mu) / f64::from(self.args.iter)

        /*
        let z_mag = (zx_sqr + zy_sqr).sqrt();
        let dz_mag = (dzx.powi(2) + dzy.powi(2)).sqrt();
        //let distance = z_mag.powi(2).ln() * z_mag / dz_mag;
        //let distance = 0. - 5. * distance.ln() / self.args.zoom.ln();

        let distance = 2. * z_mag * z_mag.ln() / dz_mag;

        //info!("distance: {}", distance);
        distance
        */
      };

      //let rgb = self.args.color_map.color(color);

      /*
      let rgb = LCH::new(
        //((color * 100. * std::f64::consts::PI).sin() / 2. + 0.5) * 100.,
        color.powf(3.5).fract() * 100.,
        0.,
        0.,
        //65.,
        //100. - (100. * color),
        //132.,//32. + 100. - (100. * color),
        //(360. * color).powi(2),
      ).rgb();
      */

      *pixel = RGB::new(
        cast(color * 255.).unwrap(),
        cast(color * 255.).unwrap(),
        cast(color * 255.).unwrap(),
      )
      .as_color();

      pp.increment();
    });

    frame
  }

  /// Creates a [`Viewport`] from [`args`](JuliaSet).
  ///
  fn viewport(args: &JuliaSet) -> Viewport {
    let w = f64::from(args.width);
    let h = f64::from(args.height);

    let aspect_ratio = w / h;

    let vp_width = aspect_ratio / args.zoom;
    let vp_height = 1. / args.zoom;

    let grid_delta_x = vp_width / w;
    let grid_delta_y = vp_height / h;

    Viewport::from_center(
      args.center.into(),
      vp_width,
      vp_height,
      grid_delta_x,
      grid_delta_y,
      args.rotation.unwrap_or(0),
    )
  }
}
