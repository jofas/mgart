use clap::{Parser, Subcommand};

use serde::Deserialize;

use std::fs::File;
use std::io::stdin;

use algorithmic_art::args::{ColorMap1dArgs, JuliaSetArgs};
use algorithmic_art::{color_map_1d, julia_set};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
  #[clap(long)]
  config: Option<String>,
  #[clap(subcommand)]
  command: Option<Command>,
}

#[derive(Subcommand, Deserialize)]
#[serde(tag = "command")]
#[serde(rename_all = "kebab-case")]
enum Command {
  JuliaSet(JuliaSetArgs),
  #[clap(name = "color-map-1d")]
  #[serde(rename = "color-map-1d")]
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

  if let Some(config) = cli.config {
    let config: Config = if config == "." {
      serde_json::from_reader(stdin()).unwrap()
    } else {
      serde_json::from_reader(File::open(config).unwrap()).unwrap()
    };

    for cmd in config.into_inner() {
      match cmd {
        Command::JuliaSet(args) => {
          println!("generating julia set with arguments:\n{}", args);
          julia_set(args);
        },
        Command::ColorMap1d(args) => {
          println!("generating 1d color map with arguments:\n{}", args);
          color_map_1d(args);
        },
      }
    }

    return;
  }

  match cli.command {
    Some(Command::JuliaSet(args)) => julia_set(args),
    Some(Command::ColorMap1d(args)) => color_map_1d(args),
    _ => println!("Successfully did nothing"),
  }
}
