use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefMutIterator,
  ParallelIterator,
};

use serde::{Deserialize, Serialize};

use display_json::DisplayAsJson;

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

  pub fn smooth(
    &self,
    buffer: &mut [f64],
    width: usize,
    height: usize,
  ) {
    let m = self.window_mean(buffer, width, height);

    let num_pixel = buffer.len();

    let pp = ProgressPrinter::new(num_pixel as u64, 100);

    for i in 0..buffer.len() {
      let x = i % width;
      let y = i / width;

      let (wx0, wy0, wx1, wy1) = discrete_rectangle_from_center(
        x as isize,
        y as isize,
        self.window_size as isize,
        self.window_size as isize,
      );

      let wx0 = wx0.max(0) as usize;
      let wy0 = wy0.max(0) as usize;
      let wx1 = wx1.min(width as isize - 1) as usize;
      let wy1 = wy1.min(height as isize - 1) as usize;

      let mut s = 0.;
      let mut cp = 0.;

      for x in wx0..=wx1 {
        for y in wy0..=wy1 {
          let j = y * width + x;

          let fpq = self.weight_function(m[i], m[j]);

          s += buffer[j] * fpq;
          cp += fpq;
        }
      }

      buffer[i] = s / cp;

      pp.increment();
    }
  }

  fn weight_function(&self, bp: f64, bq: f64) -> f64 {
    (-(bq - bp).powi(2) / self.h.powi(2)).exp()
  }

  fn window_mean(
    &self,
    buffer: &[f64],
    width: usize,
    height: usize,
  ) -> Vec<f64> {
    let sat = SummedAreaTable::new(buffer, width, height);

    let mut res = vec![0.; width * height];

    res.par_iter_mut().enumerate().for_each(|(i, p)| {
      let x = i % width;
      let y = i / width;

      *p = sat.mean_rectangle_from_center(x, y, self.n, self.n);
    });

    res
  }
}

struct SummedAreaTable {
  sat: Vec<f64>,
  width: usize,
  height: usize,
}

impl SummedAreaTable {
  fn new(buffer: &[f64], width: usize, height: usize) -> Self {
    let mut sat = buffer.to_vec();

    for i in 1..sat.len() {
      let x = i % width;
      let y = i / width;

      if x > 0 {
        sat[i] += sat[i - 1];
      }

      if y > 0 {
        sat[i] += sat[(y - 1) * width + x];
      }

      if x > 0 && y > 0 {
        sat[i] -= sat[(y - 1) * width + x - 1];
      }
    }

    Self { sat, width, height }
  }

  fn sum_rectangle_from_center(
    &self,
    cx: usize,
    cy: usize,
    size_x: usize,
    size_y: usize,
  ) -> f64 {
    let (x0, y0, x1, y1) = discrete_rectangle_from_center(
      cx as isize,
      cy as isize,
      size_x as isize,
      size_y as isize,
    );

    let x0 = x0 - 1;
    let y0 = y0 - 1;

    let x1 = x1.min(self.width as isize - 1) as usize;
    let y1 = y1.min(self.height as isize - 1) as usize;

    let c00 = if x0 >= 0 && y0 >= 0 {
      self.sat[y0 as usize * self.width + x0 as usize]
    } else {
      0.
    };

    let c01 = if y0 >= 0 {
      self.sat[y0 as usize * self.width + x1]
    } else {
      0.
    };

    let c10 = if x0 >= 0 {
      self.sat[y1 * self.width + x0 as usize]
    } else {
      0.
    };

    let c11 = self.sat[y1 * self.width + x1];

    c00 + c11 - c01 - c10
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
      cx as isize,
      cy as isize,
      size_x as isize,
      size_y as isize,
    );

    let x0 = x0.max(0);
    let y0 = y0.max(0);

    let x1 = x1.min((self.width - 1) as isize);
    let y1 = y1.min((self.height - 1) as isize);

    c / ((x1 - x0 + 1) * (y1 - y0 + 1)) as f64
  }
}

fn discrete_rectangle_from_center(
  cx: isize,
  cy: isize,
  size_x: isize,
  size_y: isize,
) -> (isize, isize, isize, isize) {
  (
    cx - size_x / 2 + isize::from(size_x % 2 == 0),
    cy - size_y / 2 + isize::from(size_y % 2 == 0),
    cx + size_x / 2,
    cy + size_y / 2,
  )
}

#[cfg(test)]
mod tests {
  use super::{
    discrete_rectangle_from_center, NonLocalMeans, SummedAreaTable,
  };

  #[test]
  fn drfc() {
    let res = discrete_rectangle_from_center(0, 0, 3, 3);
    assert_eq!(res, (-1, -1, 1, 1));

    let res = discrete_rectangle_from_center(1, 1, 3, 3);
    assert_eq!(res, (0, 0, 2, 2));

    let res = discrete_rectangle_from_center(2, 2, 3, 3);
    assert_eq!(res, (1, 1, 3, 3));

    let res = discrete_rectangle_from_center(3, 3, 3, 3);
    assert_eq!(res, (2, 2, 4, 4));

    let res = discrete_rectangle_from_center(0, 0, 2, 2);
    assert_eq!(res, (0, 0, 1, 1));

    let res = discrete_rectangle_from_center(1, 1, 2, 2);
    assert_eq!(res, (1, 1, 2, 2));

    let res = discrete_rectangle_from_center(2, 2, 2, 2);
    assert_eq!(res, (2, 2, 3, 3));

    let res = discrete_rectangle_from_center(3, 3, 2, 2);
    assert_eq!(res, (3, 3, 4, 4));
  }

  #[test]
  fn sat() {
    let image = [1.; 16];
    let width = 4;
    let height = 4;

    let sat = SummedAreaTable::new(&image, width, height);

    assert_eq!(
      sat.sat,
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

    let image: Vec<f64> = (1..=16).map(f64::from).collect();
    let width = 4;
    let height = 4;

    let sat = SummedAreaTable::new(&image, width, height);

    assert_eq!(
      sat.sat,
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
    let image: Vec<f64> = (1..=16).map(f64::from).collect();
    let width = 4;
    let height = 4;

    let sat = SummedAreaTable::new(&image, width, height);

    let res: Vec<f64> = (0..16)
      .map(|i| {
        let x = i % width;
        let y = i / width;

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
    let image: Vec<f64> = (1..=16).map(f64::from).collect();
    let width = 4;
    let height = 4;

    let sat = SummedAreaTable::new(&image, width, height);

    let res: Vec<f64> = (0..16)
      .map(|i| {
        let x = i % width;
        let y = i / width;

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
    let image: Vec<f64> = (1..=16).map(f64::from).collect();
    let width = 4;
    let height = 4;

    let nlm = NonLocalMeans::new(3, 4, 3.);

    let res = nlm.window_mean(&image, width, height);

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
  fn nlm_weight_function() {
    let nlm = NonLocalMeans::new(3, 4, 3.);

    for x in [0., 0.25, 0.33, 0.5, 0.66, 0.75, 0.99, 1.] {
      assert_eq!(nlm.weight_function(x, x), 1.);
    }
  }

  #[test]
  fn nlm() {
    let mut image: Vec<f64> = (1..=16).map(f64::from).collect();
    let width = 4;
    let height = 4;

    let nlm = NonLocalMeans::new(3, 4, 3.);

    nlm.smooth(&mut image, width, height);

    assert_eq!(
      image.into_iter().map(|x| x as u8).collect::<Vec<u8>>(),
      vec![
        [3, 4, 5, 5],
        [5, 6, 6, 7],
        [10, 11, 12, 13],
        [12, 13, 13, 14],
      ]
      .into_iter()
      .flatten()
      .collect::<Vec<u8>>(),
    );
  }
}
