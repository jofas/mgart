use std::error::Error;
use std::fmt;

/// This is thrown by operations that are disabled given the current
/// state of the system.
///
#[derive(Debug)]
pub struct OperationDisabled;

impl fmt::Display for OperationDisabled {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    write!(
      f,
      "the current state of the system does not allow this operation"
    )
  }
}

impl Error for OperationDisabled {}
