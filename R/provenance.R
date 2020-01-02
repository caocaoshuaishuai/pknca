#' Add a hash and associated information to enable checking object 
#' provenance.
#' 
#' @param object The object to add provenance
#' @param replace Replace provenance if the object already has a 
#'   provenance attribute.  (If the object already has provenance and
#'   \code{replace} is \code{FALSE}, then an error will be raised.)
#' @return The object with provenance as an added item
#' @seealso \code{\link{checkProvenance}}
#' @export
#' @importFrom digest sha1
#' @importFrom utils sessionInfo
addProvenance <- function(object, replace=FALSE) {
  if (replace) {
    attr(object, "provenance") <- NULL
  }
  if (is.null(attr(object, "provenance", exact=TRUE))) {
    # Get most of the provenance added
    prov <-
      list(
        sessionInfo=utils::sessionInfo(),
        datetime=Sys.time(),
        sysInfo=Sys.info()
      )
    class(prov) <- c("provenance", class(prov))
    prov$hash <- digest::sha1(x=c(digest::sha1(object), digest::sha1(prov)))
    attr(object, "provenance") <- prov
  } else {
    stop("object already has provenance and the option to replace it was not selected.")
  }
  object
}

#' Check the hash of an object to confirm its provenance.
#'
#' @param object The object to check provenance for
#' @return \code{TRUE} if the provenance is confirmed to be consistent, 
#'   \code{FALSE} if the provenance is not consistent, or \code{NA} if
#'   provenance is not present.
#' @seealso \code{\link{addProvenance}}
#' @export
checkProvenance <- function(object) {
  orig_provenance <- attr(object, "provenance", exact=TRUE)
  if (is.null(orig_provenance)) {
    NA
  } else {
    hash <- orig_provenance$hash
    attr(object, "provenance") <- NULL
    new_provenance <- digest::sha1(x=c(digest::sha1(object), digest::sha1(orig_provenance)))
    ret <- (hash == new_provenance)
    if (!ret) {
      message("Provenance hash mismatch: old: ", hash, ", vs. current: ", new_provenance)
    }
    ret
  }
}

#' Print the summary of a provenance object
#' 
#' @param x The object to be printed
#' @param ... Ignored
#' @return invisible text of the printed information
#' @export
print.provenance <- function(x, ...) {
  ret <- sprintf("Provenance hash %s generated on %s with %s.",
                 x$hash,
                 x$datetime,
                 x$sessionInfo$R.version[["version.string"]])
  cat(ret, "\n", sep="")
  invisible(ret)
}

#' @importFrom digest sha1
sha1.provenance <- function(x, digits=14L, zapsmall=7L, ..., algo="sha1") {
  digest::sha1(
    c(
      digest::sha1(x$sessionInfo, digits=digits, zapsmall=zapsmall, ..., algo=algo),
      digest::sha1(x$datetime, digits=digits, zapsmall=zapsmall, ..., algo=algo),
      digest::sha1(x$sysInfo, digits=digits, zapsmall=zapsmall, ..., algo=algo)
    ),
    digits=digits, zapsmall=zapsmall, ..., algo=algo
  )
}
