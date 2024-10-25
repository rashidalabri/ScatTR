use std::path::PathBuf;

use clap::{Parser, Subcommand};

use crate::commands::{
    AlignCommandArgs, DefineCommandArgs, ExtractCommandArgs, GenotypeCommandArgs, StatsCommandArgs,
};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
#[command(arg_required_else_help(true))]
pub struct Cli {
    #[arg(short, long, default_value = "1")]
    /// Number of threads to use. Only used for `depth` and `define` commands.
    pub threads: usize,

    #[arg(short, long, action = clap::ArgAction::Count)]
    pub verbosity: u8,

    /// Prefix for output files (including directory)
    #[arg()]
    pub output_prefix: PathBuf,

    #[command(subcommand)]
    pub command: Option<CliCommands>,
}

#[derive(Subcommand)]
pub enum CliCommands {
    Stats(StatsCommandArgs),
    Extract(ExtractCommandArgs),
    Define(DefineCommandArgs),
    Genotype(GenotypeCommandArgs),
    Align(AlignCommandArgs),
}
