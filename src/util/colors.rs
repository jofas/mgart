use serde::{Serialize, Deserialize};

#[derive(Clone, Copy, Serialize, Deserialize, PartialEq, Debug)]
#[serde(from = "u32")]
#[serde(into = "u32")]
pub struct Rgba {
  r: u8,
  g: u8,
  b: u8,
  a: u8,
}

impl Rgba {
  /// Creates a new [Rgba] from hex.
  ///
  /// **Example:**
  ///
  /// ```rust
  /// use algorithmic_art::util::colors::Rgba;
  ///
  /// // sunshine yellow
  /// let c = Rgba::new_hex(0xFFFD37FF);
  ///
  /// assert_eq!(c.r(), 255);
  /// assert_eq!(c.g(), 253);
  /// assert_eq!(c.b(), 55);
  /// assert_eq!(c.a(), 255);
  /// ```
  ///
  pub fn new_hex(color: u32) -> Self {
    Self {
      r: ((color & 0xFF000000) >> 24) as u8,
      g: ((color & 0xFF0000) >> 16) as u8,
      b: ((color & 0xFF00) >> 8) as u8,
      a: (color & 0xFF) as u8,
    }
  }

  /// Creates a new [Rgba] from rgba.
  ///
  /// **Example:**
  ///
  /// ```rust
  /// use algorithmic_art::util::Rgba;
  ///
  /// // sunshine yellow
  /// let c = Rgba::new_rgba(255, 253, 55, 255);
  ///
  /// assert_eq!(c.r(), 255);
  /// assert_eq!(c.g(), 253);
  /// assert_eq!(c.b(), 55);
  /// assert_eq!(c.a(), 255);
  ///
  /// assert_eq!(c.as_u32(), 0xFFFD37FF);
  /// ```
  ///
  pub fn new_rgba(r: u8, g: u8, b: u8, a: u8) -> Self {
    Self { r, g, b, a }
  }

  pub fn r(&self) -> u8 {
    self.r
  }

  pub fn g(&self) -> u8 {
    self.g
  }

  pub fn b(&self) -> u8 {
    self.b
  }

  pub fn a(&self) -> u8 {
    self.a
  }

  pub fn as_vec(&self) -> [u8; 4] {
    [self.r, self.g, self.b, self.a]
  }

  pub fn as_u32(&self) -> u32 {
    ((self.r() as u32) << 24)
      | ((self.g() as u32) << 16)
      | ((self.b() as u32) << 8)
      | (self.a() as u32)
  }
}

impl From<u32> for Rgba {
  fn from(x: u32) -> Self {
    Self::new_hex(x)
  }
}

impl Into<u32> for Rgba {
  fn into(self) -> u32 {
    self.as_u32()
  }
}

