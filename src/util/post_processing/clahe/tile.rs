use num::cast;

pub trait HistogramEqualization {
    fn transform(&self, v: f64) -> f64;
}

pub struct Tile {
    hist: Vec<u32>,
    cdf_min: u32,
    n: u32,
}

impl Tile {
    pub fn new(buffer: impl Iterator<Item = f64>, bin_count: u32, contrast_limit: u32) -> Self {
        let mut hist = vec![0; bin_count.try_into().unwrap()];
        let mut n = 0;

        // contrast limiting
        let mut clv: u32 = 0;

        for v in buffer {
            let bin = cast::<_, usize>(v * f64::from(bin_count - 1)).unwrap();

            if hist[bin] < contrast_limit {
                hist[bin] += 1;
            } else {
                clv += 1;
            }

            n += 1;
        }

        loop {
            let clv_bin: u32 = clv / bin_count;

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
        let bin = v * cast::<_, f64>(self.hist.len() - 1).unwrap();
        let bin = cast::<_, usize>(bin).unwrap();

        let res = f64::from(self.hist[bin].saturating_sub(self.cdf_min));

        res / f64::from(self.n - self.cdf_min)
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
