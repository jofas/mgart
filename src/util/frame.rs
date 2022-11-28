use map_macro::vec_no_clone;

use std::ops::{Index, IndexMut};

use crate::util::coloring::colors::Color;

/// A [`Frame`] is a two dimensional array.
///
/// A [`Frame`] is contiguous in memory with a row-first layout and
/// ungrowable.
///
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
  pub fn new(buf: Vec<T>, width: usize, height: usize) -> Self {
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

  #[must_use]
  pub fn get(&self, x: usize, y: usize) -> Option<&T> {
    if x < self.width && y < self.height {
      self.buf.get(y * self.width + x)
    } else {
      None
    }
  }

  pub fn get_mut(&mut self, x: usize, y: usize) -> Option<&mut T> {
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
}

impl<T: Default> Frame<T> {
  /// Creates a new [`Frame`] with `width` and `height`, where each
  /// element is the [Default] value of the given type.
  ///
  #[must_use]
  pub fn filled_default(width: usize, height: usize) -> Self {
    Self {
      width,
      height,
      buf: vec_no_clone![T::default(); width * height],
    }
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

    pixels.chunks_exact_mut(3).zip(self.buf.iter()).for_each(
      |(pixel, v)| {
        let rgb = v.rgb();

        pixel[0] = rgb.r();
        pixel[1] = rgb.g();
        pixel[2] = rgb.b();
      },
    );

    image::save_buffer(
      filename,
      &pixels,
      self.width as u32,
      self.height as u32,
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
