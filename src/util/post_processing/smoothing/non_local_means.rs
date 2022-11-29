use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

use divrem::DivRem;

use crate::util::frame::Frame;
use crate::util::ProgressPrinter;

#[derive(
  Serialize, Deserialize, DisplayAsJson, Clone, PartialEq, Debug,
)]
pub struct NonLocalMeans {
  n: usize,
  window_size: usize,
  h: f64,
}

impl NonLocalMeans {
  #[must_use]
  pub fn new(n: usize, window_size: usize, h: f64) -> Self {
    Self { n, window_size, h }
  }

  pub fn smooth(&self, frame: &mut Frame<f64>) {
    let width = frame.width();
    let height = frame.height();

    let m = self.window_mean(frame);

    let pp = ProgressPrinter::new(frame.len() as u64, 100);

    for i in 0..frame.len() {
      let (y, x) = i.div_rem(width);

      let (wx0, wy0, wx1, wy1) = discrete_rectangle_from_center(
        x,
        y,
        self.window_size,
        self.window_size,
        width,
        height,
      );

      let mut s = 0.;
      let mut cp = 0.;

      for x in wx0..wx1 {
        for y in wy0..wy1 {
          let fpq = self.weight_function(m[i], m[(x, y)]);

          s += frame[(x, y)] * fpq;
          cp += fpq;
        }
      }

      frame[i] = s / cp;

      pp.increment();
    }
  }

  fn weight_function(&self, bp: f64, bq: f64) -> f64 {
    (-(bq - bp).powi(2) / self.h.powi(2)).exp()
  }

  fn window_mean(&self, frame: &Frame<f64>) -> Frame<f64> {
    let sat = SummedAreaTable::new(frame);

    let mut res = Frame::filled(0., frame.width(), frame.height());

    res.par_for_each_mut(|(i, p)| {
      let (y, x) = i.div_rem(frame.width());

      *p = sat.mean_rectangle_from_center(x, y, self.n, self.n);
    });

    res
  }
}

struct SummedAreaTable(Frame<f64>);

impl SummedAreaTable {
  fn new(frame: &Frame<f64>) -> Self {
    let mut sat = frame.clone();

    for i in 1..sat.len() {
      let (y, x) = i.div_rem(frame.width());

      if x > 0 {
        sat[i] += sat[i - 1];
      }

      if y > 0 {
        sat[i] += sat[(x, y - 1)];
      }

      if x > 0 && y > 0 {
        sat[i] -= sat[(x - 1, y - 1)];
      }
    }

    Self(sat)
  }

  fn sum_rectangle_from_center(
    &self,
    cx: usize,
    cy: usize,
    size_x: usize,
    size_y: usize,
  ) -> f64 {
    let (x0, y0, x1, y1) = discrete_rectangle_from_center(
      cx,
      cy,
      size_x,
      size_y,
      self.0.width() - 1,
      self.0.height() - 1,
    );

    let mut c00 = 0.;
    let mut c01 = 0.;
    let mut c10 = 0.;

    if x0 > 0 && y0 > 0 {
      c00 = self.0[(x0 - 1, y0 - 1)];
    }

    if y0 > 0 {
      c01 = self.0[(x1, y0 - 1)];
    }

    if x0 > 0 {
      c10 = self.0[(x0 - 1, y1)];
    }

    self.0[(x1, y1)] + c00 - c01 - c10
  }

  fn mean_rectangle_from_center(
    &self,
    cx: usize,
    cy: usize,
    size_x: usize,
    size_y: usize,
  ) -> f64 {
    let c = self.sum_rectangle_from_center(cx, cy, size_x, size_y);

    let (x0, y0, x1, y1) = discrete_rectangle_from_center(
      cx,
      cy,
      size_x,
      size_y,
      self.0.width() - 1,
      self.0.height() - 1,
    );

    c / ((x1 - x0 + 1) * (y1 - y0 + 1)) as f64
  }
}

fn discrete_rectangle_from_center(
  cx: usize,
  cy: usize,
  size_x: usize,
  size_y: usize,
  upper_boundary_x: usize,
  upper_boundary_y: usize,
) -> (usize, usize, usize, usize) {
  (
    (cx + usize::from(size_x % 2 == 0)).saturating_sub(size_x / 2),
    (cy + usize::from(size_y % 2 == 0)).saturating_sub(size_y / 2),
    cx.saturating_add(size_x / 2).min(upper_boundary_x),
    cy.saturating_add(size_y / 2).min(upper_boundary_y),
  )
}

#[cfg(test)]
mod tests {
  use divrem::DivRem;

  use crate::util::frame::Frame;

  use super::{
    discrete_rectangle_from_center, NonLocalMeans, SummedAreaTable,
  };

  #[test]
  fn drfc() {
    let res = discrete_rectangle_from_center(0, 0, 3, 3, 3, 3);
    assert_eq!(res, (0, 0, 1, 1));

    let res = discrete_rectangle_from_center(1, 1, 3, 3, 3, 3);
    assert_eq!(res, (0, 0, 2, 2));

    let res = discrete_rectangle_from_center(2, 2, 3, 3, 3, 3);
    assert_eq!(res, (1, 1, 3, 3));

    let res = discrete_rectangle_from_center(3, 3, 3, 3, 3, 3);
    assert_eq!(res, (2, 2, 3, 3));

    let res = discrete_rectangle_from_center(0, 0, 2, 2, 3, 3);
    assert_eq!(res, (0, 0, 1, 1));

    let res = discrete_rectangle_from_center(1, 1, 2, 2, 3, 3);
    assert_eq!(res, (1, 1, 2, 2));

    let res = discrete_rectangle_from_center(2, 2, 2, 2, 3, 3);
    assert_eq!(res, (2, 2, 3, 3));

    let res = discrete_rectangle_from_center(3, 3, 2, 2, 3, 3);
    assert_eq!(res, (3, 3, 3, 3));
  }

  #[test]
  fn sat() {
    let frame = Frame::new(vec![1.; 16], 4, 4);

    let sat = SummedAreaTable::new(&frame);

    assert_eq!(
      sat.0.inner().to_vec(),
      vec![
        [1., 2., 3., 4.],
        [2., 4., 6., 8.],
        [3., 6., 9., 12.],
        [4., 8., 12., 16.]
      ]
      .into_iter()
      .flatten()
      .collect::<Vec<f64>>(),
    );

    let frame = Frame::new((1..=16).map(f64::from).collect(), 4, 4);

    let sat = SummedAreaTable::new(&frame);

    assert_eq!(
      sat.0.inner().to_vec(),
      vec![
        [1., 3., 6., 10.],
        [6., 14., 24., 36.],
        [15., 33., 54., 78.],
        [28., 60., 96., 136.]
      ]
      .into_iter()
      .flatten()
      .collect::<Vec<f64>>(),
    );
  }

  #[test]
  fn sat_sum_rectangle_from_center() {
    let frame = Frame::new((1..=16).map(f64::from).collect(), 4, 4);

    let sat = SummedAreaTable::new(&frame);

    let res: Vec<f64> = (0..16)
      .map(|i| {
        let (y, x) = i.div_rem(4_usize);

        sat.sum_rectangle_from_center(x, y, 3, 3)
      })
      .collect();

    assert_eq!(
      res,
      vec![
        [14., 24., 30., 22.],
        [33., 54., 63., 45.],
        [57., 90., 99., 69.],
        [46., 72., 78., 54.]
      ]
      .into_iter()
      .flatten()
      .collect::<Vec<f64>>(),
    );
  }

  #[test]
  fn sat_mean_rectangle_from_center() {
    let frame = Frame::new((1..=16).map(f64::from).collect(), 4, 4);

    let sat = SummedAreaTable::new(&frame);

    let res: Vec<f64> = (0..16)
      .map(|i| {
        let (y, x) = i.div_rem(4_usize);

        sat.mean_rectangle_from_center(x, y, 3, 3)
      })
      .collect();

    assert_eq!(
      res,
      vec![
        [3.5, 4., 5., 5.5],
        [5.5, 6., 7., 7.5],
        [9.5, 10., 11., 11.5],
        [11.5, 12., 13., 13.5]
      ]
      .into_iter()
      .flatten()
      .collect::<Vec<f64>>(),
    );
  }

  #[test]
  fn nlm_window_mean() {
    let frame = Frame::new((1..=16).map(f64::from).collect(), 4, 4);

    let nlm = NonLocalMeans::new(3, 4, 3.);

    let res = nlm.window_mean(&frame);

    assert_eq!(
      res.inner().to_vec(),
      vec![
        [3.5, 4., 5., 5.5],
        [5.5, 6., 7., 7.5],
        [9.5, 10., 11., 11.5],
        [11.5, 12., 13., 13.5]
      ]
      .into_iter()
      .flatten()
      .collect::<Vec<f64>>(),
    );
  }

  #[test]
  fn nlm_weight_function() {
    let nlm = NonLocalMeans::new(3, 4, 3.);

    for x in [0., 0.25, 0.33, 0.5, 0.66, 0.75, 0.99, 1.] {
      assert_eq!(nlm.weight_function(x, x), 1.);
    }
  }

  #[test]
  fn nlm() {
    let mut frame =
      Frame::new((1..=16).map(f64::from).collect(), 4, 4);

    let nlm = NonLocalMeans::new(3, 4, 3.);

    nlm.smooth(&mut frame);

    assert_eq!(
      frame.map(|x| x as u8).inner().to_vec(),
      vec![
        [2, 3, 4, 5],
        [4, 5, 6, 7],
        [10, 11, 12, 13],
        [12, 12, 13, 14],
      ]
      .into_iter()
      .flatten()
      .collect::<Vec<u8>>(),
    );
  }
}
