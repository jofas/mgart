use serde::{Deserialize, Serialize};

use crate::util::errors::OperationDisabled;
use crate::util::frame::Frame;

mod pos;
mod strided;
mod tile;

use pos::Pos;
use strided::Strided;
use tile::{HistogramEqualization, Tile};

#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Debug)]
pub struct CLAHE {
  contrast_limit: usize,
  bin_count: usize,
  tile_size_x: usize,
  tile_size_y: usize,
}

impl CLAHE {
  #[must_use]
  pub fn new(
    contrast_limit: usize,
    bin_count: usize,
    tile_size_x: usize,
    tile_size_y: usize,
  ) -> Self {
    Self {
      contrast_limit,
      bin_count,
      tile_size_x,
      tile_size_y,
    }
  }

  #[must_use]
  pub fn tile_size_x(&self) -> usize {
    self.tile_size_x
  }

  #[must_use]
  pub fn tile_size_y(&self) -> usize {
    self.tile_size_y
  }

  /// Applies `self` to the provided `buffer`.
  ///
  /// `buffer` has a width of `width` and height of `height`.
  ///
  /// # Panics
  ///
  /// Panics, if `width` isn't divisible by [`tile_size_x`] or
  /// `height` isn't divisible by [`tile_size_y`].
  ///
  pub fn apply(
    &self,
    buffer: &mut [f64],
    width: usize,
    height: usize,
  ) {
    // TODO: anyhow
    assert!(
      width % self.tile_size_x == 0 && height % self.tile_size_y == 0,
      "width and height must be divisible by tile_size_x and tile_size_y, respectively",
    );

    let tiles = self.tiles(buffer, width, height).unwrap();

    // TODO: parallelize?
    for (i, v) in buffer.iter_mut().enumerate() {
      let x = i % width;
      let y = i / width;

      let pos = Pos::new(
        x % self.tile_size_x,
        y % self.tile_size_y,
        self.tile_size_x,
        self.tile_size_y,
      );

      let (x_tile, y_tile) = self.tile_indices(x, y);

      if let Pos::Center = pos {
        *v = tiles.get(x_tile, y_tile).transform(*v);
        continue;
      }

      // TODO: code below into InterpolationTiles struct

      let (nw, ne, se, sw) =
        Self::interpolation_tiles(&pos, &tiles, x_tile, y_tile);

      let (dn, ds, dw, de) = self.interpolation_distances(x, y, &pos);

      let (dn, ds, dw, de) = Self::handle_corners_and_borders(
        nw, ne, se, sw, dn, ds, dw, de,
      );

      let q_nw = nw.transform(*v) * dn * dw;
      let q_ne = ne.transform(*v) * dn * de;
      let q_se = se.transform(*v) * ds * de;
      let q_sw = sw.transform(*v) * ds * dw;

      *v = q_nw + q_ne + q_se + q_sw;
    }
  }

  fn tiles(
    &self,
    buffer: &[f64],
    width: usize,
    height: usize,
  ) -> Result<Frame<Tile>, OperationDisabled> {
    let tiles_w = width / self.tile_size_x;
    let tiles_h = height / self.tile_size_y;

    let mut tiles: Frame<Tile> = Frame::empty(tiles_w, tiles_h);

    for start_block in 0..tiles_h {
      for offset in 0..tiles_w {
        let stride = Strided::new(
          width,
          self.tile_size_x,
          offset * self.tile_size_x,
          Some(self.tile_size_y),
          Some(start_block * self.tile_size_y),
          buffer,
        );

        tiles.push(Tile::new(
          stride,
          self.bin_count,
          self.contrast_limit,
        ))?;
      }
    }

    Ok(tiles)
  }

  fn interpolation_tiles<'a>(
    pos: &Pos,
    tiles: &'a Frame<Tile>,
    x: usize,
    y: usize,
  ) -> (
    Option<&'a Tile>,
    Option<&'a Tile>,
    Option<&'a Tile>,
    Option<&'a Tile>,
  ) {
    let (x1, x2, y1, y2) = match pos {
      Pos::NW => {
        (x.checked_sub(1), Some(x), y.checked_sub(1), Some(y))
      }
      Pos::NE => {
        (Some(x), x.checked_add(1), y.checked_sub(1), Some(y))
      }
      Pos::SE => {
        (Some(x), x.checked_add(1), Some(y), y.checked_add(1))
      }
      Pos::SW => {
        (x.checked_sub(1), Some(x), Some(y), y.checked_add(1))
      }
      // center pixels are handled, before this method is called
      Pos::Center => unreachable!(),
    };

    let mut nw = None;
    let mut ne = None;
    let mut se = None;
    let mut sw = None;

    if let (Some(x1), Some(y1)) = (x1, y1) {
      nw = tiles.get(x1, y1);
    }

    if let (Some(x2), Some(y1)) = (x2, y1) {
      ne = tiles.get(x2, y1);
    }

    if let (Some(x2), Some(y2)) = (x2, y2) {
      se = tiles.get(x2, y2);
    }

    if let (Some(x1), Some(y2)) = (x1, y2) {
      sw = tiles.get(x1, y2);
    }

    (nw, ne, se, sw)
  }

  fn interpolation_distances(
    &self,
    x: usize,
    y: usize,
    pos: &Pos,
  ) -> (f64, f64, f64, f64) {
    let center_x = self.tile_size_x as f64 / 2.;
    let center_y = self.tile_size_y as f64 / 2.;

    let dx = (center_x - (x % self.tile_size_x) as f64).abs();
    let dy = (center_y - (y % self.tile_size_y) as f64).abs();

    let dx = dx / self.tile_size_x as f64;
    let dy = dy / self.tile_size_y as f64;

    match pos {
      Pos::NW => (dy, 1. - dy, dx, 1. - dx),
      Pos::NE => (dy, 1. - dy, 1. - dx, dx),
      Pos::SE => (1. - dy, dy, 1. - dx, dx),
      Pos::SW => (1. - dy, dy, dx, 1. - dx),
      // center pixels are handled, before this method is called
      Pos::Center => unreachable!(),
    }
  }

  #[allow(clippy::too_many_arguments)]
  fn handle_corners_and_borders(
    nw: Option<&Tile>,
    ne: Option<&Tile>,
    se: Option<&Tile>,
    sw: Option<&Tile>,
    dn: f64,
    ds: f64,
    dw: f64,
    de: f64,
  ) -> (f64, f64, f64, f64) {
    match (nw, ne, se, sw) {
      // corners
      (Some(_), None, None, None) => (1., 0., 1., 0.),
      (None, Some(_), None, None) => (1., 0., 0., 1.),
      (None, None, Some(_), None) => (0., 1., 0., 1.),
      (None, None, None, Some(_)) => (0., 1., 1., 0.),
      // borders
      (Some(_), Some(_), None, None) => (1., 0., dw, de),
      (Some(_), None, None, Some(_)) => (dn, ds, 1., 0.),
      (None, Some(_), Some(_), None) => (dn, ds, 0., 1.),
      (None, None, Some(_), Some(_)) => (0., 1., dw, de),
      // center
      (Some(_), Some(_), Some(_), Some(_)) => (dn, ds, dw, de),
      _ => unreachable!(),
    }
  }

  fn tile_indices(&self, x: usize, y: usize) -> (usize, usize) {
    (x / self.tile_size_x, y / self.tile_size_y)
  }
}

#[cfg(test)]
mod tests {
  use super::CLAHE;

  #[test]
  fn tile_indexing() {
    let c = CLAHE::new(0, 0, 8, 8);

    assert_eq!((0, 0), c.tile_indices(0, 0));
    assert_eq!((0, 0), c.tile_indices(7, 0));
    assert_eq!((0, 0), c.tile_indices(0, 7));
    assert_eq!((0, 0), c.tile_indices(7, 7));

    assert_eq!((1, 1), c.tile_indices(8, 8));
    assert_eq!((1, 1), c.tile_indices(15, 8));
    assert_eq!((1, 1), c.tile_indices(8, 15));
    assert_eq!((1, 1), c.tile_indices(15, 15));

    assert_eq!((2, 2), c.tile_indices(16, 16));
    assert_eq!((2, 2), c.tile_indices(23, 16));
    assert_eq!((2, 2), c.tile_indices(16, 23));
    assert_eq!((2, 2), c.tile_indices(23, 23));

    assert_eq!((3, 3), c.tile_indices(24, 24));
    assert_eq!((3, 3), c.tile_indices(31, 24));
    assert_eq!((3, 3), c.tile_indices(24, 31));
    assert_eq!((3, 3), c.tile_indices(31, 31));
  }
}
