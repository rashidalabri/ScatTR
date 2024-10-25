use std::sync::Mutex;

use anyhow::{anyhow, Result};
use clap::Parser;
use log::*;
use scattr::cli::{Cli, CliCommands};
use scattr::commands::{
    run_align_command, run_define_command, run_extract_command, run_genotype_command,
    run_stats_command,
};
use slog::{o, Drain};
use slog_scope::{self, GlobalLoggerGuard};

fn main() {
    let mut cli: Cli = Cli::parse();
    let _guard = setup_logging(cli.verbosity).unwrap();
    let result = run(&mut cli);
    match result {
        Ok(_) => {
            info!("Finished");
            std::process::exit(exitcode::OK);
        }
        Err(e) => {
            error!("{}", e);
            std::process::exit(exitcode::SOFTWARE);
        }
    }
}

fn run(cli: &mut Cli) -> Result<()> {
    init_thread_pool(cli.threads)?;
    initial_debug_messages(cli.verbosity, cli.threads);
    match &cli.command {
        Some(CliCommands::Stats(args)) => run_stats_command(args, &cli.output_prefix),
        Some(CliCommands::Extract(args)) => run_extract_command(args, &cli.output_prefix),
        Some(CliCommands::Define(args)) => run_define_command(args, &cli.output_prefix),
        Some(CliCommands::Genotype(args)) => run_genotype_command(args, &cli.output_prefix),
        Some(CliCommands::Align(args)) => run_align_command(args, &cli.output_prefix),
        None => Err(anyhow!("No command was specified")),
    }?;
    Ok(())
}

fn setup_logging(verbosity: u8) -> Result<GlobalLoggerGuard> {
    let decorator = slog_term::TermDecorator::new().build();
    let drain = Mutex::new(
        slog_term::FullFormat::new(decorator)
            .use_original_order()
            .build()
            .filter_level(match verbosity {
                0 => slog::Level::Info,
                1 => slog::Level::Debug,
                2 => slog::Level::Trace,
                _ => slog::Level::Trace,
            })
            .fuse(),
    )
    .fuse();
    let logger = slog::Logger::root(drain, o!());
    let guard = slog_scope::set_global_logger(logger);
    slog_stdlog::init()?;
    Ok(guard)
}

fn init_thread_pool(threads: usize) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;
    Ok(())
}

fn initial_debug_messages(verbosity: u8, threads: usize) {
    debug!(
        "Running {} v{}",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION")
    );
    debug!("Logging verbosity set to {}", verbosity);
    debug!("Using {} threads", threads);
}
