use clap::Parser;

use serde::Deserialize;

use std::fs::File;
use std::io::stdin;

use mgart::args::{ColorMap1dArgs, JuliaSetArgs};
use mgart::buddhabrot::{buddhabrot, Args as BuddhabrotArgs};
use mgart::{color_map_1d, julia_set_interior_distance};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
  #[clap(long)]
  config: String,
}

#[derive(Deserialize)]
#[serde(tag = "command")]
#[serde(rename_all = "snake_case")]
enum Command {
  JuliaSet(JuliaSetArgs),
  Buddhabrot(BuddhabrotArgs),
  #[serde(rename = "color_map_1d")]
  ColorMap1d(ColorMap1dArgs),
}

#[derive(Deserialize)]
struct Config(Vec<Command>);

impl Config {
  fn into_inner(self) -> Vec<Command> {
    self.0
  }
}

fn main() {
  let cli = Cli::parse();

  let config: Config = if cli.config == "." {
    serde_json::from_reader(stdin()).unwrap()
  } else {
    serde_json::from_reader(File::open(cli.config).unwrap()).unwrap()
  };

  for cmd in config.into_inner() {
    match cmd {
      Command::JuliaSet(args) => {
        println!("generating julia set with arguments:\n{}", args);
        julia_set_interior_distance(args);
      }
      Command::Buddhabrot(args) => {
        println!("generating buddhabrot with arguments: \n{}", args);
        buddhabrot(args);
      }
      Command::ColorMap1d(args) => {
        println!("generating 1d color map with arguments:\n{}", args);
        color_map_1d(args);
      }
    }
  }
}
