#[derive(Debug)]
pub struct Viewport {
  pub x_min: f64,
  pub y_min: f64,
  pub x_max: f64,
  pub y_max: f64,
}

impl Viewport {
  pub fn new(x_min: f64, y_min: f64, x_max: f64, y_max: f64) -> Self {
    Self {
      x_min,
      y_min,
      x_max,
      y_max,
    }
  }

  pub fn contains_point(&self, x: f64, y: f64) -> bool {
    self.x_min <= x
      && x < self.x_max
      && self.y_min <= y
      && y < self.y_max
  }

  pub fn grid_pos(
    &self,
    x: f64,
    y: f64,
    delta_x: f64,
    delta_y: f64,
  ) -> Option<(usize, usize)> {
    if self.contains_point(x, y) {
      let x = ((x - self.x_min) / delta_x) as usize;
      let y = ((y - self.y_min) / delta_y) as usize;

      Some((x, y))
    } else {
      None
    }
  }
}

#[cfg(test)]
mod tests {
  use super::Viewport;

  #[test]
  fn grid_pos1() {
    let vp = Viewport::new(0., 0., 2., 2.);

    let (x, y) = vp.grid_pos(1.5, 0.5, 1., 1.).unwrap();

    assert_eq!(x, 1);
    assert_eq!(y, 0);
  }

  #[test]
  fn grid_pos2() {
    let vp = Viewport::new(0., 0., 2., 2.);

    let p = vp.grid_pos(1.5, -0.5, 1., 1.);

    assert_eq!(p, None);
  }

  #[test]
  fn grid_pos3() {
    let vp = Viewport::new(0., 0., 2., 2.);

    let p = vp.grid_pos(1.5, 2.5, 1., 1.);

    assert_eq!(p, None);
  }

  #[test]
  fn grid_pos4() {
    let vp = Viewport::new(-50., -50., 50., 50.);

    let (x, y) = vp.grid_pos(0., 0., 1., 1.).unwrap();

    assert_eq!(x, 50);
    assert_eq!(y, 50);
  }

  #[test]
  fn grid_pos5() {
    let vp = Viewport::new(-1., -1., 1., 1.);

    let (x, y) = vp.grid_pos(0., 0., 0.01, 0.005).unwrap();

    assert_eq!(x, 100);
    assert_eq!(y, 200);
  }

  #[test]
  fn grid_pos6() {
    let vp = Viewport::new(-1., -1., 1., 1.);

    let (x, y) = vp.grid_pos(0.999, 0.999, 0.01, 0.005).unwrap();

    assert_eq!(x, 199);
    assert_eq!(y, 399);
  }
}
