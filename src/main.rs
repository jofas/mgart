use clap::Parser;

use std::fs::File;
use std::io::stdin;

use mgart::Commands;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
  #[clap(long)]
  config: String,
}

fn main() {
  let cli = Cli::parse();

  let cmds: Commands = if cli.config == "." {
    serde_json::from_reader(stdin()).unwrap()
  } else {
    serde_json::from_reader(File::open(cli.config).unwrap()).unwrap()
  };

  cmds.execute();
}
