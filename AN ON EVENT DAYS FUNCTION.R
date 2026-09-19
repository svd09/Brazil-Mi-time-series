af_an_from_beta_event_days <- function(beta_l, x, y,
                            trunc_excess = TRUE,
                            use_mu = FALSE, mu = NULL,
                            event_mask = NULL   # optional logical: which days are “event days”
){
  n <- length(x)
  L <- length(beta_l) - 1
  if (!is.null(mu) && !use_mu) warning("mu provided but use_mu=FALSE; mu will be ignored.")
  if (use_mu && is.null(mu)) stop("Set use_mu=FALSE or supply mu (fitted means).")
  
  # --- linear predictor via lag convolution ---
  eta <- numeric(n)
  for (l in 0:L) {
    idx <- (l + 1):n
    eta[idx] <- eta[idx] + beta_l[l + 1] * x[seq_len(n - l)]
  }
  
  # --- relative risk & attributable numbers per day ---
  RR <- pmax(exp(eta), 1e-12)
  if (trunc_excess) {
    AN <- y * pmax(0, 1 - 1 / RR)    # “excess-only” attribution
  } else {
    AN <- y * (1 - 1 / RR)           # allow negative (protective) days
  }
  
  # --- totals across all days (same as your current output) ---
  AN_tot <- sum(AN, na.rm = TRUE)
  denom_all <- if (use_mu) sum(mu, na.rm = TRUE) else sum(y, na.rm = TRUE)
  AF_all <- AN_tot / denom_all
  
  # --- define event days ---
  if (is.null(event_mask)) {
    # If x looks binary, treat x==1 as event day; else require user to pass event_mask
    ux <- unique(na.omit(x))
    if (length(setdiff(ux, c(0, 1))) == 0) {
      event_mask <- (x == 1)
    } else {
      stop("Provide event_mask for non-binary x (e.g., event_mask = x >= threshold).")
    }
  } else {
    if (!is.logical(event_mask) || length(event_mask) != n)
      stop("event_mask must be a logical vector of length(x).")
  }
  
  # --- metrics restricted to event days ---
  idx_ev <- which(event_mask)
  n_event_days <- length(idx_ev)
  
  if (n_event_days > 0) {
    AN_sum_event_days   <- sum(AN[idx_ev], na.rm = TRUE)
    AN_mean_per_event   <- mean(AN[idx_ev], na.rm = TRUE)  # “AN per event-day”
    denom_event <- if (use_mu) sum(mu[idx_ev], na.rm = TRUE) else sum(y[idx_ev], na.rm = TRUE)
    AF_on_event_days    <- AN_sum_event_days / denom_event
  } else {
    AN_sum_event_days <- 0
    AN_mean_per_event <- NA_real_
    AF_on_event_days  <- NA_real_
  }
  
  # --- return a tidy summary + the daily AN series if you want to inspect ---
  out <- data.frame(
    AF_all            = AF_all,
    AN_tot_all        = AN_tot,
    n_days_all        = n,
    n_event_days      = n_event_days,
    AN_sum_event_days = AN_sum_event_days,
    AN_per_event_day  = AN_mean_per_event,
    AF_on_event_days  = AF_on_event_days
  )
  
  # attr(out, "AN_daily")      <- AN
  # attr(out, "RR_daily")      <- RR
  # attr(out, "event_mask")    <- event_mask
  return(out)
}


#- example
# beta_l = B*b
# af_an_from_beta_event_days(beta = beta_l,x = x,y=y)







