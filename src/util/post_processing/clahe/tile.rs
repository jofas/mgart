pub trait HistogramEqualization {
  fn transform(&self, v: f64) -> f64;
}

pub struct Tile {
  hist: Vec<usize>,
  cdf_min: usize,
  n: usize,
}

impl Tile {
  pub fn new(
    buffer: impl Iterator<Item = f64>,
    bin_count: usize,
    contrast_limit: usize,
  ) -> Self {
    let mut hist = vec![0; bin_count];
    let mut n = 0;

    // contrast limiting
    let mut clv = 0;

    for v in buffer {
      let bin = (v * (bin_count - 1) as f64) as usize;

      if hist[bin] < contrast_limit {
        hist[bin] += 1;
      } else {
        clv += 1;
      }

      n += 1;
    }

    loop {
      let clv_bin = clv / bin_count;

      let mut new_clv = 0;

      for b in &mut hist {
        *b += clv_bin;

        if *b > contrast_limit {
          new_clv += *b - contrast_limit;
          *b = contrast_limit;
        }
      }

      if new_clv == 0 || new_clv == clv {
        break;
      }

      clv = new_clv;
    }

    // sum up histograms
    for i in 1..hist.len() {
      hist[i] += hist[i - 1];
    }

    // minimum value in histogram
    let mut cdf_min = 0;
    for b in &hist {
      if *b > 0 {
        cdf_min = *b;
        break;
      }
    }

    Self { hist, cdf_min, n }
  }
}

impl HistogramEqualization for Tile {
  fn transform(&self, v: f64) -> f64 {
    let bin = (v * (self.hist.len() - 1) as f64) as usize;

    let res =
      self.hist[bin].checked_sub(self.cdf_min).unwrap_or(0) as f64;

    res / (self.n - self.cdf_min) as f64
  }
}

impl HistogramEqualization for Option<&Tile> {
  fn transform(&self, v: f64) -> f64 {
    match self {
      Some(t) => t.transform(v),
      None => 0.,
    }
  }
}
