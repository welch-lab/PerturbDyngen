#' Read a config
#'
#' A config file is an R script that defines a list named `cfg`.
#'
#' @param config path of a config file, the name of a bundled config (see
#'   [bundled_configs()]), or a list (returned unchanged).
#' @return the config list.
#' @export
read_config <- function(config) {
  if (is.list(config)) return(config)
  path <- config_path(config)
  e <- new.env()
  sys.source(path, envir = e)
  if (!exists("cfg", envir = e, inherits = FALSE)) stop("config must define a list named `cfg`: ", path)
  get("cfg", envir = e)
}

#' Bundled config files
#'
#' @param type "datasets" (screens of the CausalODE benchmark) or "networks" (how
#'   each bundled network was built).
#' @return named vector of file paths.
#' @export
bundled_configs <- function(type = c("datasets", "networks")) {
  type <- match.arg(type)
  f <- list.files(system.file("configs", type, package = "perturbdyngen"), pattern = "\\.R$",
                  full.names = TRUE)
  stats::setNames(f, sub("\\.R$", "", basename(f)))
}

#' Path of a config
#'
#' @param x a file path or the name of a bundled config.
#' @export
config_path <- function(x) {
  if (file.exists(x)) return(normalizePath(x))
  all <- c(bundled_configs("datasets"), bundled_configs("networks"))
  if (!x %in% names(all)) stop("no config file or bundled config named '", x, "'")
  all[[x]]
}

#' Run a config
#'
#' Builds a network (`type = "network"`) with [build_network()] or simulates a screen
#' (`type = "dataset"`) with [simulate_screen()].
#'
#' @param config path, bundled config name or list.
#' @param out_dir output directory.
#' @param ... further arguments to [simulate_screen()] or [build_network()].
#' @param overrides named list of config fields that replace those in the file.
#' @export
run_config <- function(config, out_dir, ..., overrides = list()) {
  cfg <- utils::modifyList(read_config(config), overrides)
  if (identical(cfg$type, "network")) build_network(cfg, out_dir, ...) else
    simulate_screen(cfg, out_dir, ...)
}
