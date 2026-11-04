#!/usr/bin/env Rscript

#' Parse `--key=value` arguments into a named list, converting numbers and logicals
parse_args <- function(args) {
  kv <- regmatches(args, regexpr("=", args), invert = TRUE)
  out <- lapply(kv, function(x) utils::type.convert(x[2], as.is = TRUE))
  stats::setNames(out, sub("^--", "", vapply(kv, `[`, "", 1)))
}

#' Run a config from the command line
#'
#' Usage: perturbdyngen.R <config> <out_dir> [--n_cores=N] [--experiments=A,B] [--field=value ...]
#' Other `--field=value` pairs override fields of the config.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  pos <- args[!startsWith(args, "--")]
  if (length(pos) != 2) stop("usage: perturbdyngen.R <config> <out_dir> [--n_cores=N] [--field=value ...]")
  opt <- parse_args(args[startsWith(args, "--")])
  n_cores <- if (is.null(opt$n_cores)) as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) else opt$n_cores
  exps <- if (is.null(opt$experiments)) NULL else strsplit(opt$experiments, ",")[[1]]
  opt$n_cores <- NULL; opt$experiments <- NULL
  cfg <- utils::modifyList(perturbdyngen::read_config(pos[1]), opt)
  if (identical(cfg$type, "network")) perturbdyngen::build_network(cfg, pos[2]) else
    perturbdyngen::simulate_screen(cfg, pos[2], n_cores = n_cores, experiments = exps)
}

main()
