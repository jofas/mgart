use std::fs::{self, File};
use std::path::Path;

use mgart::Algorithms;

#[test]
fn parse_examples() {
  let mut files = Vec::new();

  find_json_files_recursively(&Path::new("examples"), &mut files);

  for file in files {
    let _: Algorithms =
      serde_json::from_reader(File::open(file).unwrap()).unwrap();
  }
}

fn find_json_files_recursively(d: &Path, res: &mut Vec<Box<Path>>) {
  if d.is_dir() {
    for entry in fs::read_dir(d).unwrap() {
      let entry = entry.unwrap();

      if entry.file_type().unwrap().is_file() {
        let extension = entry.file_name().into_string().unwrap();

        let extension = extension.split(".").last().unwrap();

        if extension == "json" {
          res.push(entry.path().into_boxed_path());
        }
      } else if entry.file_type().unwrap().is_dir() {
        find_json_files_recursively(entry.path().as_path(), res);
      }
    }
  }
}
