use crate::util::errors::OperationDisabled;

/// A [`Frame`] is a two dimensional array.
///
/// A [`Frame`] is contiguous in memory with a row-first layout.
/// It is ungrowable, though it can be instantiated without all values
/// present.
/// Values can be added later using the [`push`] method, until the
/// [`Frame`] has reached its full capacity.
///
pub struct Frame<T> {
  width: usize,
  height: usize,
  buf: Vec<T>,
}

impl<T> Frame<T> {
  pub fn empty(width: usize, height: usize) -> Self {
    Self {
      width,
      height,
      buf: Vec::with_capacity(width * height),
    }
  }

  /// The width or amount of elements per row of `self`.
  ///
  pub fn width(&self) -> usize {
    self.width
  }

  /// The height or amount of rows of `self`.
  ///
  pub fn height(&self) -> usize {
    self.height
  }

  /// Returns a reference to the data buffer of `self`.
  ///
  pub fn inner(&self) -> &[T] {
    &self.buf
  }

  /// Returns a mutable reference to the data buffer of `self`.
  ///
  pub fn inner_mut(&mut self) -> &mut [T] {
    &mut self.buf
  }

  pub fn get(&self, x: usize, y: usize) -> Option<&T> {
    self.buf.get(y * self.width + x)
  }

  pub fn get_mut(&mut self, x: usize, y: usize) -> Option<&mut T> {
    self.buf.get_mut(y * self.width + x)
  }

  /// Appends a new value `v` to the end of `self`.
  ///
  /// # Errors
  ///
  /// A [`Frame`] has a capacity of [`width`] times [`height`] it can
  /// not outgrow.
  /// If you try calling this method on a [`Frame`] that is already
  /// fully allocated, [`OperationDisabled`] will be thrown.
  ///
  pub fn push(&mut self, v: T) -> Result<(), OperationDisabled> {
    if self.buf.len() == self.buf.capacity() {
      Err(OperationDisabled)
    } else {
      self.buf.push(v);
      Ok(())
    }
  }
}
