#' @keywords internal
#' @importFrom dplyr bind_rows filter distinct anti_join select
#' @importFrom tibble tibble
#' @importFrom igraph graph_from_data_frame
#' @importFrom methods as
#' @import Matrix
"_PACKAGE"

utils::globalVariables(c("from", "to", "sign"))
