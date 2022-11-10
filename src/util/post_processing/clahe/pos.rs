use std::cmp::Ordering;

#[derive(Debug)]
pub enum Pos {
  NW,
  NE,
  SE,
  SW,
  Center,
}

impl Pos {
  pub fn new(x: usize, y: usize, x_max: usize, y_max: usize) -> Self {
    let pos_x = PosH::new(x, x_max);
    let pos_y = PosV::new(y, y_max);

    match (pos_y, pos_x) {
      (PosV::N, PosH::W) => Self::NW,
      (PosV::N, PosH::E) => Self::NE,
      (PosV::N, PosH::Center) => Self::NW,
      (PosV::S, PosH::E) => Self::SE,
      (PosV::S, PosH::W) => Self::SW,
      (PosV::S, PosH::Center) => Self::SE,
      (PosV::Center, PosH::E) => Self::NE,
      (PosV::Center, PosH::W) => Self::SW,
      (PosV::Center, PosH::Center) => Self::Center,
    }
  }
}

enum PosH {
  W,
  E,
  Center,
}

impl PosH {
  fn new(x: usize, x_max: usize) -> Self {
    if x_max % 2 == 0 {
      if x < x_max / 2 {
        Self::W
      } else {
        Self::E
      }
    } else {
      match x.cmp(&(x_max / 2)) {
        Ordering::Less => Self::W,
        Ordering::Equal => Self::Center,
        Ordering::Greater => Self::E,
      }
    }
  }
}

enum PosV {
  N,
  S,
  Center,
}

impl PosV {
  fn new(y: usize, y_max: usize) -> Self {
    if y_max % 2 == 0 {
      if y < y_max / 2 {
        Self::N
      } else {
        Self::S
      }
    } else {
      match y.cmp(&(y_max / 2)) {
        Ordering::Less => Self::N,
        Ordering::Equal => Self::Center,
        Ordering::Greater => Self::S,
      }
    }
  }
}
