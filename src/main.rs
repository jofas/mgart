use clap::Parser;

use jsonnet::JsonnetVm;

use std::fs::File;
use std::io::stdin;

use mgart::Algorithms;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
  #[clap(default_value = ".")]
  file: String,
}

fn main() {
  let cli = Cli::parse();

  let cmds: Algorithms = if cli.file == "." {
    serde_json::from_reader(stdin()).unwrap()
  } else if let Some(extension) = cli.file.split(".").last() {
    match extension {
      "json" => serde_json::from_reader(File::open(cli.file).unwrap()).unwrap(),
      "jsonnet" => {
        let mut vm = JsonnetVm::new();
        let content = vm.evaluate_file(cli.file).unwrap();
        serde_json::from_str(content.as_str()).unwrap()
      }
      _ => panic!("unrecognizable file extension: {}. Only .json and .jsonnet files are supported", extension),
    }
  } else {
    panic!("unrecognizable file type. Only .json and .jsonnet files are supported");
  };

  cmds.execute();
}
