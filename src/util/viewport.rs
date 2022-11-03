use num_complex::Complex64;

use std::f64::consts::PI;

#[derive(Debug)]
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
  pub fn from_center(
    center: Complex64,
    width: f64,
    height: f64,
    grid_delta_x: f64,
    grid_delta_y: f64,
    rotation: usize,
  ) -> Self {
    let x_min = center.re - width * 0.5;
    let x_max = center.re + width * 0.5;

    let y_min = center.im - height * 0.5;
    let y_max = center.im + height * 0.5;

    let rotation = Complex64::new(
      (rotation as f64 * PI / 180.).cos(),
      (rotation as f64 * PI / 180.).sin(),
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

  pub fn contains_point(&self, p: &Complex64) -> bool {
    let p = p * self.rotation;

    self.x_min <= p.re
      && p.re < self.x_max
      && self.y_min <= p.im
      && p.im < self.y_max
  }

  pub fn grid_pos(&self, p: &Complex64) -> Option<(usize, usize)> {
    if self.contains_point(p) {
      let p = p * self.rotation;

      let x = ((p.re - self.x_min) / self.grid_delta_x) as usize;
      let y = ((p.im - self.y_min) / self.grid_delta_y) as usize;

      Some((x, y))
    } else {
      None
    }
  }
}

#[cfg(test)]
mod tests {
  use num_complex::Complex64;

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
  fn grid_pos7() {
    let center = Complex64::new(0., 0.);
    let width = 1.;
    let height = 2.;
    let grid_delta_x = 0.2;
    let grid_delta_y = 0.2;
    let rotation = 90;

    let vp = Viewport::from_center(
      center,
      width,
      height,
      grid_delta_x,
      grid_delta_y,
      rotation,
    );

    let (x, y) = vp.grid_pos(&Complex64::new(0.1, 0.)).unwrap();

    assert_eq!(x, 2);
    assert_eq!(y, 5);
  }
}
