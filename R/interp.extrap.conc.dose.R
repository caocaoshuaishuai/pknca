# Choices for events in interp.extrap.conc.dose.  This is included here to assist with testing later.
event_choices_interp.extrap.conc.dose <-
  list(conc_dose_iv_bolus_after="conc_dose_iv_bolus_after",
       conc_dose="conc_dose",
       dose_iv_bolus_after="dose_iv_bolus_after",
       dose="dose",
       conc="conc",
       output_only="output_only",
       none="none")

#' @importFrom dplyr case_when
#' @describeIn interp.extrap.conc Interpolate and extrapolate 
#'   concentrations without interpolating or extrapolating beyond doses.
#' @export
interp.extrap.conc.dose <- function(conc, time,
                                    time.dose, route.dose="extravascular", duration.dose=NA,
                                    time.out, out.after=FALSE,
                                    options=list(),
                                    conc.blq=PKNCA.choose.option("conc.blq", options),
                                    conc.na=PKNCA.choose.option("conc.na", options),
                                    ...,
                                    check=TRUE) {
  if (check) {
    check.conc.time(conc, time)
    data_conc <-
      clean.conc.blq(conc, time,
                     conc.blq=conc.blq, conc.na=conc.na,
                     check=FALSE)
  } else {
    data_conc <- data.frame(conc, time)
  }
  # Check other inputs
  if (!is.character(route.dose)) {
    route.dose <- as.character(route.dose)
  }
  if (!(all(route.dose %in% c("extravascular", "intravascular")))) {
    stop("route.dose must be either 'extravascular' or 'intravascular'")
  }
  if (!(length(route.dose) %in% c(1, length(time.dose)))) {
    stop("route.dose must either be a scalar or the same length as time.dose")
  }
  if (!all(is.na(duration.dose) | (is.numeric(duration.dose) & !is.factor(duration.dose)))) {
    stop("duration.dose must be NA or a number.")
  }
  if (!(length(duration.dose) %in% c(1, length(time.dose)))) {
    stop("duration.dose must either be a scalar or the same length as time.dose")
  }
  
  # Generate a single timeline
  
  # Concentrations are assumed to occur before dosing
  data_conc$out_after <- FALSE
  data_dose <-
    merge(
      data.frame(dose=TRUE,
                 time=time.dose,
                 route=route.dose,
                 duration=duration.dose,
                 iv_bolus=route.dose %in% "intravascular" & duration.dose %in% 0,
                 stringsAsFactors=FALSE),
      # Expand IV bolus dosing to have a before and after concentration
      data.frame(iv_bolus=c(FALSE, TRUE, TRUE),
                 out_after=c(FALSE, FALSE, TRUE)),
      all.x=TRUE)
  data_out <-
    data.frame(out=TRUE,
               out_after=out.after,
               out_order=1:length(time.out),
               time=time.out)
  data_all <-
    merge(merge(data_conc,
                data_dose,
                all=TRUE),
          data_out,
          all=TRUE)
  data_all$dose_event <- data_all$dose %in% TRUE
  data_all$conc_event <- !is.na(data_all$conc)
  data_all$iv_bolus <- data_all$iv_bolus %in% TRUE
  data_all$out <- data_all$out %in% TRUE
  data_all$dose_count <- cumsum(data_all$dose_event)
  mask_include_before <- data_all$dose_event & data_all$conc_event & !data_all$out_after
  data_all$dose_count_prev <- data_all$dose_count - mask_include_before
  data_all$event <- dplyr::case_when(
    data_all$dose_event & data_all$conc_event & data_all$iv_bolus & data_all$out_after~event_choices_interp.extrap.conc.dose$conc_dose_iv_bolus_after,
    data_all$dose_event & data_all$conc_event~event_choices_interp.extrap.conc.dose$conc_dose,
    data_all$dose_event & data_all$iv_bolus & data_all$out_after~event_choices_interp.extrap.conc.dose$dose_iv_bolus_after,
    data_all$dose_event~event_choices_interp.extrap.conc.dose$dose,
    data_all$conc_event~event_choices_interp.extrap.conc.dose$conc,
    data_all$out~event_choices_interp.extrap.conc.dose$output_only, # interpolation/extrapolation-only row
    TRUE~"unknown") # should never happen
  if (any(mask_unknown <- data_all$event %in% "unknown")) {
    # All events should be accounted for
    stop("Unknown event in interp.extrap.conc.dose at time(s): ",
         paste(unique(data_all$time[mask_unknown]), collapse=", ")) # nocov
  }
  # Remove "output_only" from event_before and event_after
  simple_locf <- function(x, missing_val) {
    mask_found <- !(x %in% missing_val)
    ret <- x[mask_found]
    ret[cumsum(mask_found)]
  }
  data_all$event_before <- simple_locf(c(event_choices_interp.extrap.conc.dose$none, data_all$event[-nrow(data_all)]),
                                       "output_only")
  data_all$event_after <- rev(simple_locf(rev(c(data_all$event[-1], event_choices_interp.extrap.conc.dose$none)),
                                          "output_only"))
  
  # Loop through the methods until all have been tested or no missing
  # values remain.
  data_all$method <- NA_character_
  for (nm in names(interp.extrap.conc.dose.select)) {
    mask <- is.na(data_all$method) &
      do.call(interp.extrap.conc.dose.select[[nm]]$select, list(x=data_all))
    if (any(mask)) {
      if ("warning" %in% names(interp.extrap.conc.dose.select[[nm]])) {
        warning(sprintf("%s: %d data points",
                        interp.extrap.conc.dose.select[[nm]]$warning,
                        sum(mask)))
        data_all$method[mask] <- nm
      } else {
        for (current_idx in which(mask)) {
          data_all$conc[current_idx] <-
            do.call(interp.extrap.conc.dose.select[[nm]]$value,
                    list(data_all=data_all,
                         current_idx=current_idx,
                         options=options,
                         ...))
          data_all$method[current_idx] <- nm
        }
      }
    }
  }
  if (any(mask_no_method <- is.na(data_all$method))) {
    # This should never happen, all eventualities should be covered
    stop("No method for imputing concentration at time(s): ",
         paste(unique(data_all$time[mask_no_method]), collapse=", ")) # nocov
  }
  # Filter to the requested time points and output
  data_out <- data_all[data_all$out,,drop=FALSE]
  data_out <- data_out[order(data_out$out_order),,drop=FALSE]
  ret <- data_out$conc
  attr(ret, "Method") <- data_out$method
  ret
}

#' @importFrom dplyr full_join bind_rows
#' @describeIn interp.extrap.conc Interpolate and extrapolate 
#'   concentrations without interpolating or extrapolating beyond doses.
#' @export
interp.extrap.conc.dose <- function(conc, time,
                                    time.dose, route.dose="extravascular", duration.dose=NA,
                                    time.out, out.after=FALSE,
                                    options=list(),
                                    conc.blq=PKNCA.choose.option("conc.blq", options),
                                    conc.na=PKNCA.choose.option("conc.na", options),
                                    ...,
                                    check=TRUE) {
  if (check) {
    check.conc.time(conc, time)
    data_conc <-
      clean.conc.blq(conc, time,
                     conc.blq=conc.blq, conc.na=conc.na,
                     check=FALSE)
  } else {
    data_conc <- data.frame(conc, time)
  }
  # Check other inputs
  if (!is.character(route.dose)) {
    route.dose <- as.character(route.dose)
  }
  if (!(all(route.dose %in% c("extravascular", "intravascular")))) {
    stop("route.dose must be either 'extravascular' or 'intravascular'")
  } else if (!(length(route.dose) %in% c(1, length(time.dose)))) {
    stop("route.dose must either be a scalar or the same length as time.dose")
  }
  if (!all(is.na(duration.dose) | (is.numeric(duration.dose) & !is.factor(duration.dose)))) {
    stop("duration.dose must be NA or a number.")
  } else if (!(length(duration.dose) %in% c(1, length(time.dose)))) {
    stop("duration.dose must either be a scalar or the same length as time.dose")
  } else if (any(na.omit(duration.dose) < 0)) {
    stop("duration.dose must be nonnegative or NA")
  }
  if (!(length(out.after) %in% c(1, length(time.out)))) {
    stop("out.after must be the same length as time.out or a scalar")
  }
  data_out <- data.frame(time=time.out,
                         after=out.after,
                         method="",
                         stringsAsFactors=FALSE)

  # Generate the timeline ####
  # In the timeline, the following columns are defined:
  # time: The current time
  # after: does the current row happen after (TRUE; as the limit from the right)
  #   or before (FALSE; as the limit from the left) the current time.  This is
  #   important for IV bolus dosing.
  # is_conc: Is the row an observed concentration?
  # is_dose: Is the row a dosing row?
  # is_dose_start: Is the row the same as a starting row for a dose
  # is_dose_ongoing: Does the row happen during an IV infusion?
  # conc: The current concentration (for concentration rows) or NA_real_
  # route: The route of dosing ("intravascular" or "extravascular")
  # duration: The duration of dosing (only permitted for route ==
  #   "intravascular")

  data_conc$is_conc <- TRUE
  data_conc$after <- FALSE
  # Generate the dosing data
  data_dose_raw <-
    data.frame(time=time.dose,
               route=route.dose,
               duration=duration.dose,
               is_dose=TRUE,
               is_dose_start=TRUE,
               after=FALSE,
               stringsAsFactors=FALSE)
  # the assumed dosing duration is 0
  data_dose_raw$duration[is.na(data_dose_raw$duration)] <- 0
  if (any(!(data_dose_raw$route %in% "intravascular") & 
          data_dose_raw$duration != 0)) {
    stop("duration.dose must be 0 except with intravascular dosing")
  }
  mask_iv <- data_dose_raw$route %in% "intravascular"
  data_dose_start <- data_dose_end <- data_dose_raw[mask_iv,]
  data_dose_both <- data_dose_raw[!mask_iv,]
  data_dose_both$after <- TRUE
  if (nrow(data_dose_start)) {
    data_dose_end$after <- TRUE
    data_dose_end$is_dose_start <- FALSE
    data_dose_end$time <- data_dose_end$time + data_dose_end$duration
  }

  # Generate the final timeline and sort it in order
  data_timeline <-
    full_join(data_conc,
              bind_rows(data_dose_start,
                        data_dose_end,
                        data_dose_both))
  # Note that using seq_len means an if block isn't needed to ensure that there
  # are any rows.
  data_timeline$is_dose_ongoing <- FALSE
  for (idx in seq_len(nrow(data_dose_start))) {
    mask_dose_ongoing <-
      (data_dose_start$time[idx] < data_timeline$time) &
      (data_timeline$time < data_dose_end$time[idx])
    data_timeline$is_dose_ongoing <-
      data_timeline$is_dose_ongoing | mask_dose_ongoing
  }
  if (any(mask_overlap_start_ongoing <-
          data_timeline$is_dose_ongoing & data_timeline$is_dose_start)) {
    stop("A dose starts while another is ongoing at time(s): ",
         paste(data_timeline$time[mask_overlap_start_ongoing], collapse=", "))
  }
  if (any(mask_overlap_end_ongoing <-
          data_timeline$is_dose_ongoing & data_timeline$after)) {
    stop("A dose ends while another is ongoing at time(s): ",
         paste(data_timeline$time[mask_overlap_end_ongoing], collapse=", "))
  }
  data_timeline <- data_timeline[order(data_timeline$time),]

  # Fill in the missing concentration data in the timeline whenever possible to
  # simplify the next step.  This is more computationally expensive than may be
  # required, but it is significantly simpler conceptually.
  data_timeline_filled <-
    iecd_predict(timeline=data_timeline, output=data_timeline, options=options, ...)
  data_out <-
    iecd_predict(timeline=data_timeline_filled, output=data_out, options=options, ...)
  data_out
}

#' @describeIn iecd_predict Find observed data
iecd_observed <- function(timeline, output, options, ...) {
  mask_observed <-
    is.na(output$conc) &
    (output$time %in% timeline$time[timeline$is_conc]) &
    (output$after %in% timeline$after[timeline$is_conc])
  for (idx in which(mask_observed)) {
    mask_idx <-
      (output$time[idx] %in% timeline$time[timeline$is_conc]) &
      (output$after[idx] %in% timeline$after[timeline$is_conc])
    output$conc[idx] <- timeline$conc[mask_idx]
  }
  output
}

#' @describeIn iecd_predict Extrapolate before the first dose of the timeline
iecd_before_beginning <- function(timeline, output, options, conc.origin=0, ...) {
  mask_before_beginning <-
    is.na(output$conc) &
    (output$time < min(timeline$time))
  if (any(mask_before_beginning)) {
    output$conc[mask_before_beginning] <- conc.origin
  }
  output
}

#' @describeIn iecd_predict Extrapolate C0
iecd_c0 <- function(timeline, output, options, ...) {
  mask_c0 <-
    is.na(output$conc) &
    (output$time %in% timeline$time[timeline$is_dose & timeline$after & timeline$duration %in% 0])
  for (idx in which(mask_c0)) {
    time_current <- output$time[idx]
    time_next_dose <- timeline$time[timeline$time > time_current &
                                      timeline$is_dose]
    if (length(time_next_dose)) {
      time_next_dose <- min(time_next_dose)
    } else {
      time_next_dose <- Inf
    }
    tmp_conc <- timeline[!is.na(timeline$conc) &
                           time_current < timeline$time &
                           (timeline$time <= time_next_dose &
                              !timeline$after),]
    if (nrow(tmp_conc)) {
      output$conc[idx] <-
        pk.calc.c0(conc=tmp_conc$conc, time=tmp_conc$time,
                   time.dose=time_current,
                   method=c("logslope", "c1"))
    } else {
      # If there is no concentration next, we cannot estimate C0.
      output$conc[idx] <- NaN
    }
  }
  output
}

#' @describeIn iecd_predict Cannot interpolate between two doses
iecd_between_doses <- function(timeline, output, options, ...) {
  mask_missing <- is.na(output$conc)
  for (idx in which(mask_missing)) {
    time_current <- output$time[idx]
    time_next_dose <- timeline$time[time_current < timeline$time &
                                      timeline$is_dose]
    time_prev_dose <- timeline$time[timeline$time < time_current &
                                      timeline$is_dose]
    if (length(time_next_dose) > 0 &
        length(time_prev_dose) > 0) {
      time_next_dose <- min(time_next_dose)
      time_prev_dose <- max(time_prev_dose)
      mask_between <-
        time_prev_dose < timeline$time &
        timeline$time < time_next_dose
      if (!any(mask_between)) {
        output$conc[idx] <- NaN
      }
    }
  }
  output
}

#' @describeIn iecd_predict Cannot extrapolate if the last event is a dose
iecd_after_dose_last <- function(timeline, output, options, ...) {
  mask_after_last <-
    is.na(output$conc) &
    timeline$is_dose[nrow(timeline)] &
    output$time > timeline$time[nrow(timeline)]
  if (any(mask_after_last)) {
    output$conc[mask_after_last] <- NaN
  }
  output
}

#' @describeIn iecd_predict Interpolate between two observations
iecd_interpolate <- function(timeline, output, options, ...) {
  mask_missing <-
    is.na(output$conc)
  for (idx in which(mask_missing)) {
    time_current <- output$time[idx]
    prior_data <- timeline[timeline$time < time_current,]
    next_data <- timeline[timeline$time > time_current,]
    if (nrow(prior_data) > 0 &
        nrow(next_data) > 0) {
      tmp_conc <-
        rbind(prior_data[nrow(prior_data),],
              next_data[1,])
      # If the prior or next data does not have a concentration, it would be NA
      # and the current concentration shouldn't be interpolated.
      output$conc[idx] <-
        interpolate.conc(conc=tmp_conc$conc, time=tmp_conc$time,
                         time.out=data_all$time[current_idx],
                         check=FALSE, ...)
    }
  }
  output
}

#' @describeIn iecd_predict Extrapolate a concentration
iecd_extrapolate <- function(timeline, output, options, ..., lambda.z) {
  mask_missing <-
    is.na(output$conc)
  for (idx in which(mask_missing)) {
    time_current <- output$time[idx]
    # Extrapolation can occur into a dose
    time_next_dose <- timeline$time[time_current <= timeline$time &
                                      timeline$is_dose]
    time_prev_dose <- timeline$time[timeline$time < time_current &
                                      timeline$is_dose]
    
  }
  output
}
  
iecd_function_order <- list(
  # Establish all methods that will not rely on other methods' results.
  "Observed concentration"=iecd_observed,
  "Before the first event (concentration or dose) in the data"=iecd_before_beginning,
  "Between two doses without a concentration"=iecd_between_doses,
  "Extrapolation after the last event in the data is a dose"=iecd_after_dose_last,
  "End of an IV bolus dose (C0)"=iecd_c0,
  # Extrapolate before interpolating so that extrapolation to the next dose time
  # can be used for interpolation.
  "Extrapolated concentration"=iecd_extrapolate,
  "Interpolated concentration"=iecd_interpolate
)

#' Do the work of interp.extrap.conc.dose.
#' @param timeline The prepared timeline (from interp.extrap.conc.dose)
#' @param output The current version of the output (from interp.extrap.conc.dose
#'   and to go back to interp.extrap.conc.dose).  Any \code{is.na(conc)} rows
#'   must be filled by the current or a subsequent method.
#' @param conc.origin The concentration to use before any observations (0 for a
#'   non-endogenous substance).
#' @param lambda.z The elimination rate for extrapolation.
#' @return \code{output} with all rows \code{!is.na(conc)}
iecd_predict <- function(timeline, output, options, ...) {
  output$method <- NA_character_
  method_idx <- 0
  while (any(is.na(output$method)) &
         method_idx < length(iecd_function_order)) {
    method_idx <- method_idx + 1
    mask_initial <- is.na(output$method)
    output <- iecd_function_order[[method_idx]](timeline, output, options, ...)
    mask_failed <- mask_initial & is.nan(output$conc)
    mask_success <- mask_initial & !mask_failed & !is.na(output$conc)
    output$method[mask_success] <- names(iecd_function_order)[method_idx]
    output$method[mask_failed] <- paste(names(iecd_function_order)[method_idx], "failed", sep=", ")
  }
  output
}

# Extrapolation ####
iecd_extrap_select <- function(x) {
  extrap_output_only <-
    x$event_before %in% c("conc_dose_iv_bolus_after", "conc") &
    x$event %in% "output_only" &
    x$event_after %in% c("dose", "none")
  extrap_dose <-
    x$event_before %in% c("conc_dose_iv_bolus_after", "conc") &
    x$event %in% "dose" &
    !(x$event_after %in% "output_only")
  extrap_output_only | extrap_dose
}
iecd_extrap_value <- function(data_all, current_idx, lambda.z, ...) {
  last_conc <- data_all[data_all$time < data_all$time[current_idx] &
                          !is.na(data_all$conc),]
  last_conc <- last_conc[nrow(last_conc),]
  if (last_conc$conc %in% 0) {
    # BLQ continues to be BLQ
    0
  } else {
    if (missing(lambda.z)) {
      lambda.z <- NA_real_
    }
    args <- list(conc=last_conc$conc[nrow(last_conc)],
                 time=last_conc$time[nrow(last_conc)],
                 time.out=data_all$time[current_idx], lambda.z=lambda.z,
                 ...)
    if (!("clast" %in% names(args))) {
      args$clast <- last_conc$conc[nrow(last_conc)]
    }
    do.call(extrapolate.conc, args)
  }
}
