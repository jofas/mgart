use serde::{Deserialize, Serialize};

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

  pub fn apply(
    &self,
    buffer: &mut [f64],
    width: usize,
    height: usize,
  ) {
    if width % self.tile_size_x != 0 || height % self.tile_size_y != 0
    {
      panic!("width and height must be divisible by tile_size_x and tile_size_y, respectively");
    }

    let tiles_w = width / self.tile_size_x;
    let tiles_h = height / self.tile_size_y;

    let tiles = self.tiles(buffer, width, tiles_w, tiles_h);

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
        let i = y_tile as usize * tiles_w + x_tile as usize;
        *v = tiles[i].transform(*v);
        continue;
      }

      let (nw, ne, se, sw) = self.interpolation_tiles(
        &pos, x_tile, y_tile, &tiles, tiles_w, tiles_h,
      );

      let (dn, ds, dw, de) = self.interpolation_distances(x, y, &pos);

      let (dn, ds, dw, de) = self.handle_corners_and_borders(
        &nw, &ne, &se, &sw, dn, ds, dw, de,
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
    tiles_w: usize,
    tiles_h: usize,
  ) -> Vec<Tile> {
    let mut tiles: Vec<Tile> = Vec::with_capacity(tiles_w * tiles_h);

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
        ));
      }
    }

    tiles
  }

  fn interpolation_tiles<'a>(
    &self,
    pos: &Pos,
    x_tile: isize,
    y_tile: isize,
    tiles: &'a [Tile],
    tiles_w: usize,
    tiles_h: usize,
  ) -> (
    Option<&'a Tile>,
    Option<&'a Tile>,
    Option<&'a Tile>,
    Option<&'a Tile>,
  ) {
    let (x1, x2, y1, y2) = match pos {
      Pos::NW => (x_tile - 1, x_tile, y_tile - 1, y_tile),
      Pos::NE => (x_tile, x_tile + 1, y_tile - 1, y_tile),
      Pos::SE => (x_tile, x_tile + 1, y_tile, y_tile + 1),
      Pos::SW => (x_tile - 1, x_tile, y_tile, y_tile + 1),
      Pos::Center => panic!("Cannot handle pixels at tile center"),
    };

    let mut nw = None;
    let mut ne = None;
    let mut se = None;
    let mut sw = None;

    if x1 >= 0 && y1 >= 0 {
      nw = Some(&tiles[y1 as usize * tiles_w + x1 as usize]);
    }

    if x2 < tiles_w as isize && y1 >= 0 {
      ne = Some(&tiles[y1 as usize * tiles_w + x2 as usize]);
    }

    if x2 < tiles_w as isize && y2 < tiles_h as isize {
      se = Some(&tiles[y2 as usize * tiles_w + x2 as usize]);
    }

    if x1 >= 0 && y2 < tiles_h as isize {
      sw = Some(&tiles[y2 as usize * tiles_w + x1 as usize]);
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
      Pos::Center => panic!("Unable to handle pixels at tile center"),
    }
  }

  #[allow(clippy::too_many_arguments)]
  fn handle_corners_and_borders(
    &self,
    nw: &Option<&Tile>,
    ne: &Option<&Tile>,
    se: &Option<&Tile>,
    sw: &Option<&Tile>,
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
      _ => panic!("impossible state"),
    }
  }

  fn tile_indices(&self, x: usize, y: usize) -> (isize, isize) {
    (
      (x / self.tile_size_x) as isize,
      (y / self.tile_size_y) as isize,
    )
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
