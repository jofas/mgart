use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use num_complex::Complex64;

use map_macro::vec_no_clone;

use log::info;

use std::sync::atomic::{AtomicU64, Ordering};

use crate::util::coloring::colors::Color;
use crate::util::coloring::ColorMap1d;
use crate::util::frame::Frame;
use crate::util::post_processing::PostProcessing;
use crate::util::sampler::{Distribution, Sampler, Sampling};
use crate::util::viewport::Viewport;
use crate::util::{ComplexNumber, ProgressPrinter};

#[derive(Serialize, Deserialize, DisplayAsJsonPretty, Clone)]
pub struct Buddhabrot {
  width: usize,
  height: usize,
  center: ComplexNumber,
  zoom: f64,
  iter: u64,
  rotation: Option<usize>,
  color_map: ColorMap1d,
  exponent: f64,
  sample_count: u64,
  sampler: Sampler,
  post_processing: Vec<PostProcessing>,
}

impl Buddhabrot {
  /// Creates a new instance of [`Buddhabrot`].
  ///
  #[must_use]
  #[allow(clippy::too_many_arguments)]
  pub fn new(
    width: usize,
    height: usize,
    center: ComplexNumber,
    zoom: f64,
    iter: u64,
    rotation: Option<usize>,
    color_map: ColorMap1d,
    exponent: f64,
    sample_count: u64,
    sampler: Sampler,
    post_processing: Vec<PostProcessing>,
  ) -> Self {
    Self {
      width,
      height,
      center,
      zoom,
      iter,
      rotation,
      color_map,
      exponent,
      sample_count,
      sampler,
      post_processing,
    }
  }

  /// Transforms `self` into a [`BuddhabrotCreator`] that can be used
  /// to create a [buddhabrot][buddhabrot] rendering.
  ///
  /// [buddhabrot]: https://en.wikipedia.org/wiki/Buddhabrot
  ///
  #[must_use]
  pub fn creator(self) -> Creator {
    Creator::new(self)
  }
}

#[derive(Clone)]
pub struct Creator {
  args: Buddhabrot,
  viewport: Viewport,
}

impl Creator {
  /// Creates a new instance of [`Creator`].
  ///
  #[must_use]
  pub fn new(args: Buddhabrot) -> Self {
    let viewport = Self::viewport(&args);

    Self { args, viewport }
  }

  /// Creates a rendering of a [buddhabrot][buddhabrot] as a `PNG`
  /// image.
  ///
  /// # Errors
  ///
  /// Returns an error, if the configuration is faulty or if the
  /// generated `PNG` image could not be saved to disk.
  ///
  /// [buddhabrot]: https://en.wikipedia.org/wiki/Buddhabrot
  ///
  #[must_use]
  pub fn create(&self) -> Frame<Color> {
    let sampler = self.sampler();

    info!("starting buddhabrot generation");

    let buf = vec_no_clone![
      AtomicU64::new(0);
      self.args.width * self.args.height
    ];

    let frame = Frame::new(buf, self.args.width, self.args.height);

    let pp = ProgressPrinter::new(self.args.sample_count, 2500);

    (0..self.args.sample_count).into_par_iter().for_each(|_| {
      let c = sampler.sample();

      let (grid_pos, j) = Self::trace_point(
        c,
        self.args.iter,
        self.args.exponent,
        &self.viewport,
      );

      if j != self.args.iter {
        for pos in grid_pos {
          frame[pos].fetch_add(1, Ordering::Relaxed);
        }
      }

      pp.increment();
    });

    info!("buddhabrot generation finished");

    info!("starting post processing");

    let mut frame = frame.map(|x| x.into_inner() as f64);

    for process in &self.args.post_processing {
      process.apply(&mut frame);
    }

    info!("post processing done");

    frame.map(|x| self.args.color_map.color(x).as_color())
  }

  fn trace_point(
    c: Complex64,
    iter: u64,
    exponent: f64,
    viewport: &Viewport,
  ) -> (Vec<(usize, usize)>, u64) {
    let mut z = c;
    let mut z_sqr = z.norm_sqr();

    let mut grid_pos = Vec::new();
    let mut j = 0;

    if let Some(pos) = viewport.rotated_grid_pos(&z) {
      grid_pos.push(pos);
    }

    while j < iter && z_sqr <= 4.0 {
      z = z.powf(exponent) + c;
      z_sqr = z.norm_sqr();

      if let Some(pos) = viewport.rotated_grid_pos(&z) {
        grid_pos.push(pos);
      }

      j += 1;
    }

    (grid_pos, j)
  }

  #[must_use]
  pub fn sampler(&self) -> Distribution<Complex64> {
    self.args.sampler.distribution(&|c| {
      let (grid_pos, iter) = Self::trace_point(
        *c,
        self.args.iter,
        self.args.exponent,
        &self.viewport,
      );

      if iter == self.args.iter {
        0.
      } else {
        grid_pos.len() as f64 / self.args.iter as f64
      }
    })
  }

  /// Creates a [`Viewport`] from [`args`](Buddhabrot).
  ///
  fn viewport(args: &Buddhabrot) -> Viewport {
    let aspect_ratio = args.width as f64 / args.height as f64;

    let vp_width = aspect_ratio / args.zoom;
    let vp_height = 1. / args.zoom;

    let grid_delta_x = vp_width / args.width as f64;
    let grid_delta_y = vp_height / args.height as f64;

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
