pub struct Strided<'a, T> {
  block_size: usize,
  elements: usize,
  offset: usize,
  block_count: Option<usize>,
  start_block: usize,
  buffer: &'a [T],
  block: usize,
  element: usize,
}

impl<'a, T> Strided<'a, T> {
  pub fn new(
    block_size: usize,
    elements: usize,
    offset: usize,
    block_count: Option<usize>,
    start_block: Option<usize>,
    buffer: &'a [T],
  ) -> Self {
    Self {
      block_size,
      elements,
      offset,
      block_count,
      start_block: start_block.unwrap_or(0),
      buffer,
      block: 0,
      element: 0,
    }
  }
}

impl<'a, T: Copy> Iterator for Strided<'a, T> {
  type Item = T;

  fn next(&mut self) -> Option<Self::Item> {
    // early exit if block_count is defined and has been reached
    match self.block_count {
      Some(bc) if bc == self.block => return None,
      _ => (),
    }

    let i = (self.start_block + self.block) * self.block_size
      + self.offset
      + self.element;

    if self.element < self.elements - 1 {
      self.element += 1;
    } else {
      self.element = 0;
      self.block += 1;
    }

    self.buffer.get(i).copied()
  }
}

#[cfg(test)]
mod tests {
  use super::Strided;

  #[test]
  fn strided() {
    let buf: Vec<f64> = (0..16).map(f64::from).collect();

    let s = Strided::new(4, 2, 0, None, None, &buf);
    let v: Vec<f64> = s.collect();
    assert_eq!(v, vec![0., 1., 4., 5., 8., 9., 12., 13.]);

    let s = Strided::new(4, 2, 2, None, None, &buf);
    let v: Vec<f64> = s.collect();
    assert_eq!(v, vec![2., 3., 6., 7., 10., 11., 14., 15.]);

    let s = Strided::new(4, 3, 1, None, None, &buf);
    let v: Vec<f64> = s.collect();
    assert_eq!(
      v,
      vec![1., 2., 3., 5., 6., 7., 9., 10., 11., 13., 14., 15.]
    );

    let s = Strided::new(4, 2, 0, Some(2), None, &buf);
    let v: Vec<f64> = s.collect();
    assert_eq!(v, vec![0., 1., 4., 5.]);

    let s = Strided::new(4, 2, 2, Some(2), None, &buf);
    let v: Vec<f64> = s.collect();
    assert_eq!(v, vec![2., 3., 6., 7.]);

    let s = Strided::new(4, 2, 0, Some(2), Some(2), &buf);
    let v: Vec<f64> = s.collect();
    assert_eq!(v, vec![8., 9., 12., 13.]);

    let s = Strided::new(4, 2, 2, Some(2), Some(2), &buf);
    let v: Vec<f64> = s.collect();
    assert_eq!(v, vec![10., 11., 14., 15.]);
  }
}
