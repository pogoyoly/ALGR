#' Assign crop by table
#'
#' @param output_obj A ALGR output object
#' @param field_crop_tbl a table connecting fields to crops
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
#'assignments <- data.frame(
#'field = c(1, 3, 7, 10),   # field indices in map$field_list
#'crop  = c(2, 1, 3, 2)     # your crop ids
#')
#'
#'output2 <- assign_crops_by_table(map, assignments)
#'return_by_crop(output2, method = 1)
assign_crops_by_table <- function(output_obj,
                                  field_crop_tbl,
                                  on_dup    = c("last","first","error"),
                                  strict    = TRUE) {

  on_dup <- match.arg(on_dup)

  field_col = NULL
  crop_col  = NULL

  # --- coerce and pick columns
  if (!is.data.frame(field_crop_tbl)) {
    stop("'field_crop_tbl' must be a data.frame-like object.")
  }
  if (is.null(field_col)) field_col <- names(field_crop_tbl)[1]
  if (is.null(crop_col))  crop_col  <- names(field_crop_tbl)[2]
  if (!(field_col %in% names(field_crop_tbl))) stop("Column '", field_col, "' not found.")
  if (!(crop_col  %in% names(field_crop_tbl)))  stop("Column '", crop_col,  "' not found.")

  df <- data.frame(field = field_crop_tbl[[field_col]],
                   crop  = field_crop_tbl[[crop_col]],
                   stringsAsFactors = FALSE)

  # --- basic checks
  n_fields <- length(output_obj$field_list)
  if (n_fields == 0) stop("output_obj$field_list appears to be empty.")

  suppressWarnings({
    df$field <- as.integer(df$field)
    df$crop  <- as.integer(df$crop)
  })

  # identify invalid rows
  bad_na   <- is.na(df$field) | is.na(df$crop)
  bad_rng  <- !(df$field >= 1 & df$field <= n_fields)
  bad_any  <- bad_na | bad_rng

  if (any(bad_any)) {
    msg <- sprintf(
      "Found %d invalid assignment row(s): %d NA field/crop, %d out-of-range field(s) (valid: 1..%d).",
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

    # keep first or last occurrence per field
    ord <- if (on_dup == "first") !duplicated(df$field) else !duplicated(df$field, fromLast = TRUE)
    df <- df[ord, , drop = FALSE]
  }

  # --- build a full crop assignment (default 0)
  crop_assign <- rep(0L, n_fields)
  crop_assign[df$field] <- df$crop

  # --- assign to all fields
  for (i in seq_len(n_fields)) {
    output_obj$field_list[[i]]@crop <- crop_assign[i]
  }

  return(output_obj)
}


