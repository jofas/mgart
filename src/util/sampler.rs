use serde::{Deserialize, Serialize};

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use num::complex::Complex64;

use rand::random;

use log::info;

use num::cast;

use std::f64::consts::PI;
use std::marker::PhantomData;

use crate::util::ProgressPrinter;

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "type")]
pub enum Sampler {
    /// Uniform random distribution generating values in the range of
    /// `[0,h]`.
    ///
    Uniform { h: f64 },
    /// Uniform random distribution generating values with a radius
    /// in `[0,r]` from the origin of a space, e.g. `(0,0)` in an
    /// unspecified 2d space or `0 + 0i` in the complex space.
    ///
    UniformPolar { r: f64 },
    /// Kernel Density Estimation.
    ///
    /// *NOTE:* `kernel` and `pre_sampler` must eventually resolve to a
    /// sampler that does not require pre-sampling, e.g.
    /// [Self::Uniform] or [Self::UniformPolar].
    ///
    KernelDensityEstimation {
        weighted: bool,
        kernel: Box<Sampler>,
        population: u64,
        p_min: f64,
        pre_sampler: Box<Sampler>,
    },
}

impl Sampler {
    /// Creates a [Distribution] from this instance of [Sampler].
    ///
    /// The `true_probability` function is needed for the
    /// [`Self::KernelDensityEstimation`] to create samples from the
    /// true probability function.
    ///
    pub fn distribution<T: Sync + Send + Copy, P: Fn(&T) -> f64 + Sync + Send>(
        &self,
        true_probability: &P,
    ) -> Distribution<T>
    where
        Distribution<T>: Sampling<Space = T>,
    {
        match self {
            Self::Uniform { h } => Distribution::Uniform(Uniform::<T>::new(*h)),
            Self::UniformPolar { r } => Distribution::UniformPolar(UniformPolar::<T>::new(*r)),
            Self::KernelDensityEstimation {
                weighted,
                kernel,
                population,
                p_min,
                pre_sampler,
            } => {
                let elems = Self::pre_sample(
                    &pre_sampler.distribution(true_probability),
                    true_probability,
                    *population,
                    *p_min,
                );

                Distribution::KernelDensityEstimation(KDE::<T>::new(
                    *weighted,
                    Box::new(kernel.distribution(true_probability)),
                    *population,
                    elems,
                ))
            }
        }
    }

    fn pre_sample<T: Sync + Send, P: Fn(&T) -> f64 + Sync + Send>(
        pre_sampler: &Distribution<T>,
        true_probability: &P,
        population: u64,
        p_min: f64,
    ) -> Vec<(T, f64)>
    where
        Distribution<T>: Sampling<Space = T>,
    {
        info!("initializing kde population");

        let pp = ProgressPrinter::new(population, 2500);

        let res = (0..population)
            .into_par_iter()
            .filter_map(|_| {
                let sample = pre_sampler.sample();

                let p = true_probability(&sample);

                pp.increment();

                if p >= p_min {
                    Some((sample, p))
                } else {
                    None
                }
            })
            .collect();

        info!("initializing kde population done");

        res
    }
}

pub trait Sampling {
    type Space;

    fn sample(&self) -> Self::Space;
    fn kernel_sample(&self, center: &Self::Space) -> Self::Space;
}

pub enum Distribution<T> {
    Uniform(Uniform<T>),
    UniformPolar(UniformPolar<T>),
    KernelDensityEstimation(KDE<T>),
}

impl Sampling for Distribution<Complex64> {
    type Space = Complex64;

    fn sample(&self) -> Self::Space {
        match self {
            Self::Uniform(u) => u.sample(),
            Self::UniformPolar(up) => up.sample(),
            Self::KernelDensityEstimation(kde) => kde.sample(),
        }
    }

    fn kernel_sample(&self, center: &Self::Space) -> Self::Space {
        match self {
            Self::Uniform(u) => u.kernel_sample(center),
            Self::UniformPolar(up) => up.kernel_sample(center),
            Self::KernelDensityEstimation(kde) => kde.kernel_sample(center),
        }
    }
}

pub struct Uniform<T> {
    h: f64,
    phantom: PhantomData<T>,
}

impl<T> Uniform<T> {
    #[must_use]
    pub fn new(h: f64) -> Self {
        Self {
            h,
            phantom: PhantomData,
        }
    }
}

impl Sampling for Uniform<Complex64> {
    type Space = Complex64;

    /// Generates a uniformly random complex number with real and
    /// imaginary parts in the `[0, h]` range.
    ///
    fn sample(&self) -> Self::Space {
        Complex64::new(random::<f64>() * self.h, random::<f64>() * self.h)
    }

    /// Generates a uniformly random complex number with real and
    /// imaginary parts in the `[center - h/2, center + h/2]` range.
    ///
    fn kernel_sample(&self, center: &Self::Space) -> Self::Space {
        let re = (random::<f64>() - 0.5) * self.h;
        let im = (random::<f64>() - 0.5) * self.h;

        Complex64::new(center.re + re, center.im + im)
    }
}

impl Sampling for Uniform<f64> {
    type Space = f64;

    /// Generates a uniformly random real number in the `[0, h]` range.
    ///
    fn sample(&self) -> Self::Space {
        random::<f64>() * self.h
    }

    /// Generates a uniformly random real number in the
    /// `[center - h/2, center + h/2]` range.
    ///
    fn kernel_sample(&self, center: &Self::Space) -> Self::Space {
        center + (random::<f64>() - 0.5) * self.h
    }
}

pub struct UniformPolar<T> {
    r: f64,
    phantom: PhantomData<T>,
}

impl<T> UniformPolar<T> {
    #[must_use]
    pub fn new(r: f64) -> Self {
        Self {
            r,
            phantom: PhantomData,
        }
    }
}

impl Sampling for UniformPolar<Complex64> {
    type Space = Complex64;

    /// Generates a uniformly random complex number with a radius `r`
    /// from the origin of the complex space `0 + 0i`.
    ///
    fn sample(&self) -> Self::Space {
        Complex64::from_polar(self.r * random::<f64>(), 2. * PI * random::<f64>())
    }

    fn kernel_sample(&self, center: &Self::Space) -> Self::Space {
        center + self.sample()
    }
}

/// Kernel Density Estimation.
///
pub struct KDE<T> {
    elems: Vec<T>,
    kernel: Box<Distribution<T>>,
}

impl<T: Copy> KDE<T> {
    /// Creates a new instance of [`KDE`].
    ///
    /// # Panics
    ///
    /// Panics, if `population` overflows the 53 bit mantissa of [f64].
    ///
    #[must_use]
    pub fn new(
        weighted: bool,
        kernel: Box<Distribution<T>>,
        population: u64,
        elems: Vec<(T, f64)>,
    ) -> Self {
        let population = cast::<_, f64>(population).unwrap();

        let elems = if weighted {
            let p_sum: f64 = elems.iter().map(|x| x.1).sum();

            let mut res = Vec::new();

            for (e, p) in elems {
                let count = cast::<_, usize>(p * population / p_sum).unwrap();

                for _ in 0..count {
                    res.push(e);
                }
            }

            res
        } else {
            let (elems, _): (Vec<T>, Vec<f64>) = elems.into_iter().unzip();
            elems
        };

        Self { elems, kernel }
    }
}

impl Sampling for KDE<Complex64> {
    type Space = Complex64;

    fn sample(&self) -> Self::Space {
        let len = cast::<_, f64>(self.elems.len()).unwrap();

        let idx = cast::<_, usize>(random::<f64>() * len).unwrap();

        self.kernel.kernel_sample(&self.elems[idx])
    }

    fn kernel_sample(&self, center: &Self::Space) -> Self::Space {
        center + self.sample()
    }
}

#[cfg(test)]
mod tests {
    use super::{Distribution, Uniform, KDE};

    #[test]
    fn weighted_kde_population() {
        let kde = KDE::new(
            true,
            Box::new(Distribution::Uniform(Uniform::new(1.))),
            10,
            vec![(0, 0.), (1, 0.25), (2, 0.25), (3, 0.5)],
        );

        assert_eq!(kde.elems, vec![1, 1, 2, 2, 3, 3, 3, 3, 3],);
    }
}
