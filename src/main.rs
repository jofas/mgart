use clap::{Parser, Subcommand};

use std::fs::File;

use algorithmic_art::args::{Config, JuliaSetArgs};
use algorithmic_art::julia_set;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
  #[clap(long)]
  config: Option<String>,
  #[clap(subcommand)]
  command: Commands,
}

#[derive(Subcommand)]
enum Commands {
  JuliaSet(JuliaSetArgs),
}

fn main() {
  let cli = Cli::parse();

  if let Some(config) = cli.config {
    let file = File::open(config).unwrap();

    let config: Config = serde_json::from_reader(file).unwrap();

    for args in config.into_inner() {
      julia_set(args);
    }

    return;
  }

  match cli.command {
    Commands::JuliaSet(args) => julia_set(args),
  }
}
