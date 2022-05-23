use clap::{Parser, Subcommand};

use std::fs::File;
use std::io::stdin;

use algorithmic_art::args::{ColorMap1dArgs, Config, JuliaSetArgs};
use algorithmic_art::{color_map_1d, julia_set};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
  #[clap(long)]
  config: Option<String>,
  #[clap(subcommand)]
  command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
  JuliaSet(JuliaSetArgs),
  #[clap(name = "color-map-1d")]
  ColorMap1d(ColorMap1dArgs),
}

fn main() {
  let cli = Cli::parse();

  if let Some(config) = cli.config {
    let config: Config = if config == "." {
      serde_json::from_reader(stdin()).unwrap()
    } else {
      serde_json::from_reader(File::open(config).unwrap()).unwrap()
    };

    for args in config.into_inner() {
      println!("generating julia set with arguments:\n{}", args);
      julia_set(args);
    }

    return;
  }

  match cli.command {
    Some(Commands::JuliaSet(args)) => julia_set(args),
    Some(Commands::ColorMap1d(args)) => color_map_1d(args),
    _ => println!("Successfully did nothing"),
  }
}
