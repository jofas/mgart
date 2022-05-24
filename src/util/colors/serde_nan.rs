use serde::de::{Deserializer, Error, Visitor};
use serde::ser::Serializer;

use std::fmt;

pub(super) fn serialize<S>(f: &f64, s: S) -> Result<S::Ok, S::Error>
where
  S: Serializer,
{
  if f.is_nan() {
    S::serialize_str(s, "NaN")
  } else {
    S::serialize_f64(s, *f)
  }
}

pub(super) fn deserialize<'de, D>(d: D) -> Result<f64, D::Error>
where
  D: Deserializer<'de>,
{
  struct StringOrNum;

  impl<'de> Visitor<'de> for StringOrNum {
    type Value = f64;

    fn expecting(
      &self,
      formatter: &mut fmt::Formatter,
    ) -> fmt::Result {
      formatter.write_str("string or number")
    }

    fn visit_str<E: Error>(self, s: &str) -> Result<f64, E> {
      if s == "NaN" {
        Ok(f64::NAN)
      } else {
        Err(Error::custom("expected a number"))
      }
    }

    fn visit_f32<E: Error>(self, f: f32) -> Result<f64, E> {
      Ok(f as f64)
    }

    fn visit_f64<E: Error>(self, f: f64) -> Result<f64, E> {
      Ok(f)
    }

    fn visit_i32<E: Error>(self, i: i32) -> Result<f64, E> {
      Ok(i as f64)
    }

    fn visit_i64<E: Error>(self, i: i64) -> Result<f64, E> {
      Ok(i as f64)
    }

    fn visit_u32<E: Error>(self, u: u32) -> Result<f64, E> {
      Ok(u as f64)
    }

    fn visit_u64<E: Error>(self, u: u64) -> Result<f64, E> {
      Ok(u as f64)
    }
  }

  d.deserialize_any(StringOrNum)
}
