use serde::Deserialize;

use log::debug;

pub mod buddhabrot;
pub mod debug;
pub mod julia_set;
pub mod util;

use crate::buddhabrot::Buddhabrot;
use crate::debug::ColorMap1dRenderer;
use crate::julia_set::JuliaSet;

#[derive(Deserialize)]
pub struct Algorithm {
    #[serde(flatten)]
    algorithm: AlgorithmInner,
    filename: String,
}

#[derive(Deserialize)]
#[serde(tag = "algorithm")]
#[serde(rename_all = "snake_case")]
pub enum AlgorithmInner {
    JuliaSet(JuliaSet),
    Buddhabrot(Buddhabrot),
    #[serde(rename = "debug.color_map_1d")]
    ColorMap1dRenderer(ColorMap1dRenderer),
}

impl Algorithm {
    /// Creates a new instance of [`Algorithm`] creating a rendering
    /// of a [`JuliaSet`].
    ///
    #[must_use]
    pub fn julia_set(julia_set: JuliaSet, filename: String) -> Self {
        Self {
            algorithm: AlgorithmInner::JuliaSet(julia_set),
            filename,
        }
    }

    /// Creates a new instance of [`Algorithm`] creating a rendering
    /// of a [`Buddhabrot`].
    ///
    #[must_use]
    pub fn buddhabrot(buddhabrot: Buddhabrot, filename: String) -> Self {
        Self {
            algorithm: AlgorithmInner::Buddhabrot(buddhabrot),
            filename,
        }
    }

    /// Creates a new instance of [`Algorithm`] creating a debug
    /// rendering of a [`ColorMap1d`][ColorMap1dRenderer].
    ///
    #[must_use]
    pub fn color_map_1d_renderer(
        color_map_1d_renderer: ColorMap1dRenderer,
        filename: String,
    ) -> Self {
        Self {
            algorithm: AlgorithmInner::ColorMap1dRenderer(color_map_1d_renderer),
            filename,
        }
    }

    /// Executes the rendering process for the given [`Algorithm`],
    /// creating a media file containing the generated artwork.
    ///
    /// # Panics
    ///
    /// Panics if the rendering process fails.
    /// Rendering processes fail, because saving the generated image to
    /// disk was unsuccessful or because the provided
    /// [`configuration`](Self) is faulty.
    ///
    pub fn create(self) {
        match self.algorithm {
            AlgorithmInner::JuliaSet(j) => {
                debug!("generating julia:\n{}", j);
                j.creator().create().save_as_image(&self.filename);
            }
            AlgorithmInner::Buddhabrot(b) => {
                debug!("generating buddhabrot: \n{}", b);
                b.creator().create().save_as_image(&self.filename);
            }
            AlgorithmInner::ColorMap1dRenderer(c) => {
                debug!("generating 1d color map:\n{}", c);
                c.creator().create().save_as_image(&self.filename);
            }
        }
    }
}

#[derive(Deserialize)]
pub struct Algorithms(Vec<Algorithm>);

impl Algorithms {
    /// Creates a new instance of [`Algorithms`].
    ///
    #[must_use]
    pub fn new(algorithms: Vec<Algorithm>) -> Self {
        Self(algorithms)
    }

    /// Executes each [`Algorithm`] successively.
    ///
    /// Multi-threading is implemented inside the rendering process of
    /// each [`Algorithm`].
    ///
    /// # Panics
    ///
    /// If one of the provided algorithms fails.
    ///
    pub fn create(self) {
        for cmd in self.0 {
            cmd.create();
        }
    }
}
