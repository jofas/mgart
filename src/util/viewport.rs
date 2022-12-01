use num::complex::Complex64;

use num::cast;

use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct Viewport {
  x_min: f64,
  y_min: f64,
  x_max: f64,
  y_max: f64,
  grid_delta_x: f64,
  grid_delta_y: f64,
  rotation: Complex64,
}

impl Viewport {
  #[must_use]
  pub fn from_center(
    center: Complex64,
    width: f64,
    height: f64,
    grid_delta_x: f64,
    grid_delta_y: f64,
    rotation: u16,
  ) -> Self {
    let x_min = center.re - width * 0.5;
    let x_max = center.re + width * 0.5;

    let y_min = center.im - height * 0.5;
    let y_max = center.im + height * 0.5;

    let rotation = Complex64::new(
      (f64::from(rotation) * PI / 180.).cos(),
      (f64::from(rotation) * PI / 180.).sin(),
    );

    Self {
      x_min,
      y_min,
      x_max,
      y_max,
      grid_delta_x,
      grid_delta_y,
      rotation,
    }
  }

  #[must_use]
  pub fn contains_point(&self, p: &Complex64) -> bool {
    self.x_min <= p.re
      && p.re < self.x_max
      && self.y_min <= p.im
      && p.im < self.y_max
  }

  #[must_use]
  pub fn contains_rotated_point(&self, p: &Complex64) -> bool {
    self.contains_point(&(p * self.rotation))
  }

  /// Takes a complex number `p` and returns the grid cell it is part
  /// of.
  ///
  /// If `p` does not lie in the [`Viewport`], [None] is returned.
  ///
  /// # Panics
  ///
  /// Could panic, if something goes wrong with casting between
  /// [f64] and [usize].
  ///
  #[must_use]
  pub fn grid_pos(&self, p: &Complex64) -> Option<(usize, usize)> {
    if self.contains_point(p) {
      let x =
        cast::<_, usize>((p.re - self.x_min) / self.grid_delta_x)
          .unwrap();

      let y =
        cast::<_, usize>((p.im - self.y_min) / self.grid_delta_y)
          .unwrap();

      Some((x, y))
    } else {
      None
    }
  }

  #[must_use]
  pub fn rotated_grid_pos(
    &self,
    p: &Complex64,
  ) -> Option<(usize, usize)> {
    self.grid_pos(&(p * self.rotation))
  }

  /// Takes a 2d grid index of `x` and `y`, returning the complex
  /// number at the grid position.
  ///
  /// # Panics
  ///
  /// Panics, if `x` or `y` overflow the 53 bit mantissa of [f64].
  ///
  #[must_use]
  pub fn point(&self, x: usize, y: usize) -> Complex64 {
    let x = cast::<_, f64>(x).unwrap();
    let y = cast::<_, f64>(y).unwrap();

    let re = x * self.grid_delta_x + self.x_min;
    let im = y * self.grid_delta_y + self.y_min;

    Complex64::new(re, im)
  }

  #[must_use]
  pub fn rotated_point(&self, x: usize, y: usize) -> Complex64 {
    self.point(x, y) * self.rotation
  }
}

#[cfg(test)]
mod tests {
  use num::complex::Complex64;

  use super::Viewport;

  #[test]
  fn grid_pos1() {
    let center = Complex64::new(1., 1.);
    let width = 2.;
    let height = 2.;
    let grid_delta_x = 1.;
    let grid_delta_y = 1.;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      0,
    );

    let (x, y) = vp.grid_pos(&Complex64::new(1.5, 0.5)).unwrap();

    assert_eq!(x, 1);
    assert_eq!(y, 0);
  }

  #[test]
  fn grid_pos2() {
    let center = Complex64::new(1., 1.);
    let width = 2.;
    let height = 2.;
    let grid_delta_x = 1.;
    let grid_delta_y = 1.;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      0,
    );

    let p = vp.grid_pos(&Complex64::new(1.5, -0.5));

    assert_eq!(p, None);
  }

  #[test]
  fn grid_pos3() {
    let center = Complex64::new(1., 1.);
    let width = 2.;
    let height = 2.;
    let grid_delta_x = 1.;
    let grid_delta_y = 1.;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      0,
    );

    let p = vp.grid_pos(&Complex64::new(1.5, 2.5));

    assert_eq!(p, None);
  }

  #[test]
  fn grid_pos4() {
    let center = Complex64::new(0., 0.);
    let width = 100.;
    let height = 100.;
    let grid_delta_x = 1.;
    let grid_delta_y = 1.;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      0,
    );

    let (x, y) = vp.grid_pos(&Complex64::new(0., 0.)).unwrap();

    assert_eq!(x, 50);
    assert_eq!(y, 50);
  }

  #[test]
  fn grid_pos5() {
    let center = Complex64::new(0., 0.);
    let width = 2.;
    let height = 2.;
    let grid_delta_x = 0.01;
    let grid_delta_y = 0.005;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      0,
    );

    let (x, y) = vp.grid_pos(&Complex64::new(0., 0.)).unwrap();

    assert_eq!(x, 100);
    assert_eq!(y, 200);
  }

  #[test]
  fn grid_pos6() {
    let center = Complex64::new(0., 0.);
    let width = 2.;
    let height = 2.;
    let grid_delta_x = 0.01;
    let grid_delta_y = 0.005;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      0,
    );

    let (x, y) = vp.grid_pos(&Complex64::new(0.999, 0.999)).unwrap();

    assert_eq!(x, 199);
    assert_eq!(y, 399);
  }

  #[test]
  fn rotated_grid_pos() {
    let center = Complex64::new(0., 0.);
    let width = 1.;
    let height = 2.;
    let grid_delta_x = 0.1;
    let grid_delta_y = 0.1;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      0,
    );

    let (x, y) =
      vp.rotated_grid_pos(&Complex64::new(0.15, 0.)).unwrap();

    assert_eq!(x, 6);
    assert_eq!(y, 10);

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      90,
    );

    let (x, y) =
      vp.rotated_grid_pos(&Complex64::new(0.15, 0.)).unwrap();

    assert_eq!(x, 5);
    assert_eq!(y, 11);

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      180,
    );

    let (x, y) =
      vp.rotated_grid_pos(&Complex64::new(0.15, 0.)).unwrap();

    assert_eq!(x, 3);
    assert_eq!(y, 10);

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      270,
    );

    let (x, y) =
      vp.rotated_grid_pos(&Complex64::new(0.15, 0.)).unwrap();

    assert_eq!(x, 5);
    assert_eq!(y, 8);
  }

  #[test]
  fn point1() {
    let center = Complex64::new(0., 0.);
    let width = 2.;
    let height = 2.;
    let grid_delta_x = 0.01;
    let grid_delta_y = 0.01;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      0,
    );

    let point = vp.point(100, 100);

    assert!(point.re <= f64::EPSILON);
    assert!(point.im <= f64::EPSILON);
  }

  #[test]
  fn point2() {
    let center = Complex64::new(0., 0.);
    let width = 2.;
    let height = 2.;
    let grid_delta_x = 0.01;
    let grid_delta_y = 0.01;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      0,
    );

    let point = vp.point(0, 0);

    assert!((point.re + 1.) <= f64::EPSILON);
    assert!((point.im + 1.) <= f64::EPSILON);
  }

  #[test]
  fn point2_with_rotation() {
    let center = Complex64::new(0., 0.);
    let width = 2.;
    let height = 2.;
    let grid_delta_x = 0.01;
    let grid_delta_y = 0.01;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      180,
    );

    let point = vp.rotated_point(0, 0);

    assert!((point.re - 1.) <= f64::EPSILON);
    assert!((point.im - 1.) <= f64::EPSILON);
  }

  #[test]
  fn point3() {
    let center = Complex64::new(0., 0.);
    let width = 2.;
    let height = 2.;
    let grid_delta_x = 0.01;
    let grid_delta_y = 0.005;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      0,
    );

    let point = vp.point(199, 399);

    assert!((point.re - 0.99).abs() <= f64::EPSILON);
    assert!((point.im - 0.995).abs() <= f64::EPSILON);
  }
}
