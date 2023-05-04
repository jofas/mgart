use serde::{Deserialize, Serialize};

use display_json::DisplayAsJsonPretty;

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use num::complex::Complex64;

use map_macro::vec_no_clone;

use log::info;

use std::sync::atomic::{AtomicU32, Ordering};

use crate::util::coloring::colors::Color;
use crate::util::coloring::ColorMap1d;
use crate::util::frame::Frame;
use crate::util::post_processing::PostProcessing;
use crate::util::sampler::{Distribution, Sampler, Sampling};
use crate::util::viewport::Viewport;
use crate::util::{ComplexNumber, ProgressPrinter};

#[derive(Serialize, Deserialize, DisplayAsJsonPretty, Clone)]
pub struct Buddhabrot {
    width: u32,
    height: u32,
    center: ComplexNumber,
    zoom: f64,
    iter: u32,
    rotation: Option<u16>,
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
        width: u32,
        height: u32,
        center: ComplexNumber,
        zoom: f64,
        iter: u32,
        rotation: Option<u16>,
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

    /// Transforms `self` into a [`Creator`] that can be used
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
    /// # Panics
    ///
    /// Panics, if the configuration is faulty, e.g if there happens an
    /// overflow.
    ///
    /// [buddhabrot]: https://en.wikipedia.org/wiki/Buddhabrot
    ///
    #[must_use]
    pub fn create(&self) -> Frame<Color> {
        let sampler = self.sampler();

        info!("starting buddhabrot generation");

        let buf = vec_no_clone![
            AtomicU32::new(0);
            (self.args.width * self.args.height).try_into().unwrap()
        ];

        let frame = Frame::new(buf, self.args.width, self.args.height);

        let pp = ProgressPrinter::new(self.args.sample_count, 2500);

        (0..self.args.sample_count).into_par_iter().for_each(|_| {
            let (grid_pos, j) = self.trace_point(sampler.sample());

            if j != self.args.iter {
                for pos in grid_pos {
                    frame[pos].fetch_add(1, Ordering::Relaxed);
                }
            }

            pp.increment();
        });

        info!("buddhabrot generation finished");

        info!("starting post processing");

        let mut frame = frame.map(|x| f64::from(x.into_inner()));

        for process in &self.args.post_processing {
            process.apply(&mut frame);
        }

        info!("post processing done");

        frame.map(|x| self.args.color_map.color(x).as_color())
    }

    #[must_use]
    pub fn sampler(&self) -> Distribution<Complex64> {
        self.args.sampler.distribution(&|c| self.probability(*c))
    }

    fn trace_point(&self, c: Complex64) -> (Vec<(usize, usize)>, u32) {
        let mut z = c;
        let mut z_sqr = z.norm_sqr();

        let mut grid_pos = Vec::with_capacity(self.args.iter as usize);
        let mut j = 0;

        if let Some(pos) = self.viewport.rotated_grid_pos(&z) {
            grid_pos.push(pos);
        }

        while j < self.args.iter && z_sqr <= 4.0 {
            z = z.powf(self.args.exponent) + c;
            z_sqr = z.norm_sqr();

            if let Some(pos) = self.viewport.rotated_grid_pos(&z) {
                grid_pos.push(pos);
            }

            j += 1;
        }

        (grid_pos, j)
    }

    fn probability(&self, c: Complex64) -> f64 {
        let mut z = c;
        let mut z_sqr = z.norm_sqr();

        let mut vp_hits = 0;
        let mut j = 0;

        if self.viewport.contains_rotated_point(&z) {
            vp_hits += 1;
        }

        while j < self.args.iter && z_sqr <= 4.0 {
            z = z.powf(self.args.exponent) + c;
            z_sqr = z.norm_sqr();

            if self.viewport.contains_rotated_point(&z) {
                vp_hits += 1;
            }

            j += 1;
        }

        if j == self.args.iter {
            0.
        } else {
            f64::from(vp_hits) / f64::from(self.args.iter)
        }
    }

    /// Creates a [`Viewport`] from [`args`](Buddhabrot).
    ///
    fn viewport(args: &Buddhabrot) -> Viewport {
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
