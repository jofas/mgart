use clap::{Parser, Subcommand};

use algorithmic_art::args::JuliaSetArgs;
use algorithmic_art::julia_set;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
  #[clap(subcommand)]
  command: Commands,
}

#[derive(Subcommand)]
enum Commands {
  JuliaSet(JuliaSetArgs),
}

fn main() {
  let cli = Cli::parse();

  match cli.command {
    Commands::JuliaSet(args) => julia_set(args),
  }
}
