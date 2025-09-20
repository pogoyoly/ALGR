#' Title
#'
#' @param output_obj A ALGR output object
#' @param field_farmer_tbl a table connecting fields to farmers
#' @param on_dup  How to handle duplicate rows for the same field with "last", "first", "error"
#' @param strict  If TRUE, stop on any invalid rows; if FALSE, drop invalid rows with a warning.
#'
#' @return A ALGR output object with crops assigned
#' @export
#'
#' @examples
#'r<-terra::rast(matrix(1, nrow=50, ncol=50))
#'map<-establish_pac(potential_space= r,
#'                   cell_size=1,
#'                   includsion_value = 1,
#'                   mean_field_size = 100,
#'                   sd_field_size = 25,
#'                   distribution = "norm",
#'                   mean_shape_index = 1,
#'                   sd_shape_index = 0.3,
#'                   percent = 75)
#'
#'set.seed(1245)
#'
#'assign_tbl <- data.frame(field = c(1, 3, 7, 10), farmer = c(2, 1, 3, 2))
#'
#'map2 <- assign_farmers_by_table(
#'  output_obj       = map,
#'  field_farmer_tbl = assign_tbl)
#'
#'
#'
#'
#'return_by_farmer(map2, method = 1)
#'
assign_farmers_by_table <- function(output_obj,
                                    field_farmer_tbl,
                                    on_dup          = c("last","first","error"),
                                    strict          = TRUE) {


  fill_unassigned = 0L
  field_col       = NULL
  farmer_col      = NULL


  on_dup <- match.arg(on_dup)

  # --- coerce and pick columns
  if (!is.data.frame(field_farmer_tbl)) {
    stop("'field_farmer_tbl' must be a data.frame-like object.")
  }
  if (is.null(field_col))  field_col  <- names(field_farmer_tbl)[1]
  if (is.null(farmer_col)) farmer_col <- names(field_farmer_tbl)[2]
  if (!(field_col %in% names(field_farmer_tbl)))  stop("Column '", field_col,  "' not found.")
  if (!(farmer_col %in% names(field_farmer_tbl))) stop("Column '", farmer_col, "' not found.")

  df <- data.frame(
    field  = field_farmer_tbl[[field_col]],
    farmer = field_farmer_tbl[[farmer_col]],
    stringsAsFactors = FALSE
  )

  # --- basic checks
  n_fields <- length(output_obj$field_list)
  if (n_fields == 0) stop("output_obj$field_list appears to be empty.")

  suppressWarnings({
    df$field  <- as.integer(df$field)
    df$farmer <- as.integer(df$farmer)
  })

  # identify invalid rows
  bad_na  <- is.na(df$field) | is.na(df$farmer)
  bad_rng <- !(df$field >= 1 & df$field <= n_fields)
  bad_any <- bad_na | bad_rng

  if (any(bad_any)) {
    msg <- sprintf(
      "Found %d invalid assignment row(s): %d NA field/farmer, %d out-of-range field(s) (valid: 1..%d).",
      sum(bad_any), sum(bad_na), sum(bad_rng), n_fields
    )
    if (strict) stop(msg)
    warning(msg, call. = FALSE)
    df <- df[!bad_any, , drop = FALSE]
  }

  # handle duplicates
  if (any(duplicated(df$field))) {
    dup_fields <- sort(unique(df$field[duplicated(df$field)]))
    msg <- sprintf("Duplicate assignments for field(s): %s.", paste(dup_fields, collapse = ", "))
    if (on_dup == "error") stop(msg)
    warning(paste(msg, "Using", on_dup, "occurrence."), call. = FALSE)

    keep <- if (on_dup == "first") !duplicated(df$field) else !duplicated(df$field, fromLast = TRUE)
    df <- df[keep, , drop = FALSE]
  }

  # --- build full assignment (default fill_unassigned) and write back
  farmer_assign <- rep(as.integer(fill_unassigned), n_fields)
  if (nrow(df) > 0) {
    farmer_assign[df$field] <- df$farmer
  }

  for (i in seq_len(n_fields)) {
    output_obj$field_list[[i]]@farmer <- farmer_assign[i]
  }

  return(output_obj)
}



