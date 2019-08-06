#' Create PKNCAconc objects in a way that is easier to program.
#'
#' @param data An object that can be coerced to a tibble with `as_tibble()`.
#' @param time_nominal,time_actual The column name (as a character string) for
#'   nominal or actual time for the concentration measurement or the start of an
#'   interval measurement.
#' @param duration The column name (as a character string) for the duration of
#'   an interval measurement.
#' @param concentration,volume,amount The column name (as a character string)
#'   for the measurement.
#' @param groups Zero or more column names for grouping variables (such as
#'   treatment)
#' @param subject The column name (as a character string) for the column
#'   indicating subject identifier.  If missing or `NULL`, data are assumed to
#'   be single-subject.
#' @param analyte The column name (as a character string) for the column
#'   indicating the analyte identifier.
#' @param exclude The column name (as a character string) for the column
#'   indicating rows of data to exclude from all calculations.
#' @param exclude_half.life,include_half.life The column name (as a character
#'   string) indicating rows to exclude from half life calculation or as the
#'   only rows to include in half-life calculations.  Both may be provided, but
#'   only one may be used for a given interval for a given subject.
#' @return A PKNCAconc object
#' @export
PKNCAconc_api <- function(data,
                          time_nominal, time_actual, duration,
                          concentration, volume, amount,
                          groups, subject, analyte,
                          exclude, exclude_half.life, include_half.life) {
  # Verify that the columns exist and set them in the column definition.
  column_definition <- list()
  for (nm in setdiff(names(formals()), "data")) {
    browser()
    if (missing(nm)) {
      column_definition[[nm]] <- NULL
    } else {
      column_definition[[nm]] <- get(x=nm)
    }
  }
  as_PKNCAconc(
    list(
      data=data,
      columns=column_definition
    )
  )
}

#' Convert an object to a PKNCAconc object after confirming that it is valid.
#'
#' @param object The object to convert or check.
#' @return The PKNCAconc object
#' @export
as_PKNCAconc <- function(object)
  UseMethod("as_PKNCAconc")

#' @importFrom tibble tibble
as_PKNCAconc.list <- function(object) {
  class(object) <- "PKNCAconc"
  if (!identical(names(object), c("data", "columns"))) {
    stop("PKNCAconc object must only contain 'data' and 'columns'.")
  }
  ## The data must have... data
  if (nrow(object$data) == 0) {
    stop("`data` must have at least one row.")
  }
  # Ensure that data is a tibble.
  object$data <- tibble::as_tibble(object$data)
  # Maximum number of columns for each definition (defaults to 1, if not
  # specified)
  max_ncols <- c("groups"=Inf)
  # Check the edits of all columns.
  allowed_classes <-
    list(
      time_nominal="numeric",
      time_actual=gsub(
        pattern="^time_calc.",
        replacement="",
        x=as.character(methods(time_calc)),
        fixed=TRUE
      ),
      duration="numeric",
      concentration="numeric",
      volume="numeric",
      amount="numeric",
      groups="*any*",
      subject="*any*",
      analyte="*any*",
      exclude="logical",
      exclude_half.life="logical",
      include_half.life="logical"
    )
  # Verify that this is in sync with PKNCAconc_api
  stopifnot(names(allowed_classes) == setdiff(names(formals(PKNCAconc_api)), "data"))
  extra_columns <-
    setdiff(names(object$columns), names(allowed_classes))
  missing_columns <-
    setdiff(names(allowed_classes), names(object$columns))
  if (length(extra_columns)) {
    warning(
      "The following column definitions are extra within the PKNCAconc object and will be ignored: ",
      paste("`", extra_columns, "`", collapse=", ")
    )
  }
  if (length(missing_columns)) {
    stop("The following column definitions are missing within the PKNCAconc object: ",
         paste("`", missing_columns, "`", collapse=", ")
    )
  }
  for (nm in names(allowed_classes)) {
    current_columns <- get_columns(object, nm)
    if (!is.null(current_columns)) {
      current_max_ncol <-
        if (nm %in% names(max_ncols)) {
          max_ncols[[nm]]
        } else {
          1
        }
      if (ncol(current_columns) > current_max_ncol) {
        stop(
          "`", nm, "` has ", ncol(current_columns),
          " columns defined, but only ", current_max_ncol,
          " columns are allowed."
        )
      }
      good_class <-
        if (any(allowed_classes[[nm]] %in% "*any*")) {
          TRUE
        } else {
          if (any(allowed_classes[[nm]] %in% "numeric") & is.numeric(current_columns[[1]])) {
            TRUE
          } else {
            any(class(current_columns[[1]]) %in% allowed_classes[[nm]])
          }
        }
      if (!good_class) {
        stop(
          "The column for ", nm, " is class ",
          paste(class(current_columns[[1]]), collapse=", "),
          " but it must be one of the following classes: ",
          paste(allowed_classes[[nm]], collapse=", ")
        )
      }
    }
  }
  if (is.null(get_columns(object, "subject"))) {
    message("subject column was not provided, assuming that the data are single-subject.")
  }
  if (is.null(get_columns(object, "concentration")) &
      is.null(get_columns(object, "volume")) &
      is.null(get_columns(object, "amount"))) {
    stop("One or more of `concentration`, `volume`, and `amount` must be given.")
  }
  if (is.null(get_columns(object, "time_nominal")) &
      is.null(get_columns(object, "time_actual"))) {
    stop("At least one time measure (`time_nominal` or `time_actual`) must be given.")
  }
}

get_columns <- function(object, name)
  UseMethod("get_columns")

get_columns.default <- function(object, name) {
  col_name <- object$columns[[name]]
  if (is.null(col_name)) {
    NULL
  } else if (!all(col_name %in% names(object$data))) {
    missing_cols <- setdiff(col_name, names(object$data))
    stop(
      "The following columns defining ", name, " are not in the data: ",
      paste0("`", missing_cols, "`", collapse=", ")
    )
  } else {
    object$data[, col_name, drop=FALSE]
  }
}

get_key_columns <- function(object)
  UseMethod("get_key_columns")

get_key_columns.PKNCAconc <- function(object, ..., time_type=c("nominal", "actual")) {
  time_type <- match.arg(time_type)
  time_col <-
    if (time_type %in% "nominal") {
      get_columns(object, "time_nominal")
    } else if (time_type %in% "actual") {
      get_columns(object, "time_actual")
    } else {
      stop("Unknown time_type") #nocov
    }
  bind_cols(
    get_columns(object, "groups"),
    get_columns(object, "subject"),
    get_columns(object, "analyte"),
    time_col
  )
}

get_key_columns.PKNCAdose <- function(object, ..., time_type=c("nominal", "actual")) {
  time_type <- match.arg(time_type)
  time_col <-
    if (time_type %in% "nominal") {
      get_columns(object, "time_nominal")
    } else if (time_type %in% "actual") {
      get_columns(object, "time_actual")
    } else {
      stop("Unknown time_type") #nocov
    }
  bind_cols(
    get_columns(object, "groups"),
    get_columns(object, "subject"),
    time_col
  )
}
