use clap::{Parser, Subcommand, Args};

use algorithmic_art::util::ColorMap1D;

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

#[derive(Args)]
struct JuliaSetArgs {
  width: u32,
  height: u32,
  zoom: f32,
  zpx: f32,
  zpy: f32,
  iter: u32,
  filename: String,
  color_map: ColorMap1D,
}

fn main() {
  let cli = Cli::parse();

  match &cli.command {
    Commands::JuliaSet(_) => println!("julia set"),
  }
}
