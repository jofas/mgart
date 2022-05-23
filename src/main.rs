use clap::{Parser, Subcommand};

use algorithmic_art::args::JuliaSetArgs;

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

  match &cli.command {
    Commands::JuliaSet(_) => println!("julia set"),
  }
}
