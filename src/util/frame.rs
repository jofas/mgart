use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

use map_macro::vec_no_clone;

use std::ops::{Index, IndexMut};

use crate::util::coloring::colors::Color;

/// A [`Frame`] is a two dimensional array.
///
/// A [`Frame`] is contiguous in memory with a row-first layout and
/// ungrowable.
///
#[derive(Debug, Clone)]
pub struct Frame<T> {
    width: usize,
    height: usize,
    buf: Vec<T>,
}

impl<T> Frame<T> {
    /// Creates a new [`Frame`] from the provided `buf` with `width` and
    /// `height`.
    ///
    /// # Panics
    ///
    /// Panics, if `buf` does not contain exactly `width` times `height`
    /// many elements.
    ///
    #[must_use]
    pub fn new<I: TryInto<usize>>(buf: Vec<T>, width: I, height: I) -> Self
    where
        <I as TryInto<usize>>::Error: std::fmt::Debug,
    {
        let width = width.try_into().unwrap();
        let height = height.try_into().unwrap();

        assert_eq!(
            buf.len(),
            width * height,
            "the provided buffer must contain width * height many elements",
        );

        Self { width, height, buf }
    }

    /// The width or amount of elements per row of `self`.
    ///
    #[must_use]
    pub fn width(&self) -> usize {
        self.width
    }

    /// The height or amount of rows of `self`.
    ///
    #[must_use]
    pub fn height(&self) -> usize {
        self.height
    }

    /// How many elements the [`Frame`] contains (`[width] * [height]`).
    ///
    #[must_use]
    pub fn len(&self) -> usize {
        self.buf.len()
    }

    /// Whether the [`Frame`] has a length of zero.
    ///
    /// ```
    /// use mgart::util::frame::Frame;
    ///
    /// assert!(Frame::<i32>::new(vec![], 0, 0).is_empty());
    /// assert!(!Frame::<i32>::new(vec![1], 1, 1).is_empty());
    /// ```
    ///
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.buf.is_empty()
    }

    /// Tries to retrieve a reference to the element at row `x` and
    /// column `y`.
    ///
    /// Returns [`None`] if either `x` or `y` are out of bounds.
    ///
    /// # Panics
    ///
    /// Panics, if either `x` or `y` do not fit into the [`usize`] type
    /// of the machine.
    ///
    #[must_use]
    pub fn get<I: TryInto<usize>>(&self, x: I, y: I) -> Option<&T>
    where
        <I as TryInto<usize>>::Error: std::fmt::Debug,
    {
        let x = x.try_into().unwrap();
        let y = y.try_into().unwrap();

        if x < self.width && y < self.height {
            self.buf.get(y * self.width + x)
        } else {
            None
        }
    }

    /// Tries to retrieve a mutable reference to the element at row `x`
    /// and column `y`.
    ///
    /// Returns [`None`] if either `x` or `y` are out of bounds.
    ///
    /// # Panics
    ///
    /// Panics, if either `x` or `y` do not fit into the [`usize`] type
    /// of the machine.
    ///
    pub fn get_mut<I: TryInto<usize>>(&mut self, x: I, y: I) -> Option<&mut T>
    where
        <I as TryInto<usize>>::Error: std::fmt::Debug,
    {
        let x = x.try_into().unwrap();
        let y = y.try_into().unwrap();

        if x < self.width && y < self.height {
            self.buf.get_mut(y * self.width + x)
        } else {
            None
        }
    }

    #[must_use]
    pub fn inner(&self) -> &[T] {
        &self.buf
    }

    pub fn inner_mut(&mut self) -> &mut [T] {
        &mut self.buf
    }

    pub fn for_each_mut(&mut self, f: impl Fn((usize, &mut T))) {
        self.buf.iter_mut().enumerate().for_each(f);
    }

    pub fn map<U>(self, f: impl Fn(T) -> U) -> Frame<U> {
        let buf = self.buf.into_iter().map(f).collect();

        Frame::new(buf, self.width, self.height)
    }
}

impl<T: Send> Frame<T> {
    pub fn par_for_each_mut<F>(&mut self, f: F)
    where
        F: Fn((usize, &mut T)) + Sync + Send,
    {
        self.buf.par_iter_mut().enumerate().for_each(f);
    }
}

impl<T: Clone> Frame<T> {
    /// Creates a new [`Frame`] with `width` and `height`, where each
    /// element is a clone of `v`.
    ///
    pub fn filled(v: T, width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            buf: vec![v; width * height],
        }
    }
}

impl<T: Default> Frame<T> {
    /// Creates a new [`Frame`] with `width` and `height`, where each
    /// element is the [Default] value of the given type.
    ///
    /// # Panics
    ///
    /// Panics, if there is an overflow (e.g. `width * height` doesn't
    /// fit in the `usize` type of the machine).
    ///
    #[must_use]
    pub fn filled_default<I: TryInto<usize>>(width: I, height: I) -> Self
    where
        <I as TryInto<usize>>::Error: std::fmt::Debug,
    {
        let width = width.try_into().unwrap();
        let height = height.try_into().unwrap();

        let buf = vec_no_clone![T::default(); width.checked_mul(height).unwrap()];

        Self { width, height, buf }
    }
}

impl Frame<Color> {
    /// Saves this [`Frame`] as an [`RGB`](RGB) image to the file with
    /// `filename`.
    ///
    /// # Panics
    ///
    /// Panics if the image file could not be created.
    ///
    /// [RGB]: crates::util::coloring::colors::RGB
    ///
    pub fn save_as_image(&self, filename: &str) {
        let mut pixels = vec![0_u8; self.len() * 3];

        pixels
            .chunks_exact_mut(3)
            .zip(self.buf.iter())
            .for_each(|(pixel, v)| {
                let rgb = v.rgb();

                pixel[0] = rgb.r();
                pixel[1] = rgb.g();
                pixel[2] = rgb.b();
            });

        image::save_buffer(
            filename,
            &pixels,
            self.width.try_into().unwrap(),
            self.height.try_into().unwrap(),
            image::ColorType::Rgb8,
        )
        .unwrap();
    }
}

impl<T> Index<usize> for Frame<T> {
    type Output = T;

    fn index(&self, i: usize) -> &Self::Output {
        &self.buf[i]
    }
}

impl<T> IndexMut<usize> for Frame<T> {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.buf[i]
    }
}

impl<T> Index<(usize, usize)> for Frame<T> {
    type Output = T;

    fn index(&self, (x, y): (usize, usize)) -> &Self::Output {
        &self.buf[y * self.width + x]
    }
}

impl<T> IndexMut<(usize, usize)> for Frame<T> {
    fn index_mut(&mut self, (x, y): (usize, usize)) -> &mut Self::Output {
        &mut self.buf[y * self.width + x]
    }
}
