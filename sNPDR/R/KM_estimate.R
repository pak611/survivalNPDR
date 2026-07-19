#=========================================================================#
#' kmImputeTimes
#'
#' Kaplan-Meier conditional mean imputation for right-censored survival times.
#'
#' For each censored subject with observed time \eqn{t_c}, replaces \eqn{t_c}
#' with the conditional mean \eqn{\hat{E}[T \mid T > t_c]}, estimated by
#' rectangular integration of the Kaplan-Meier step function:
#' \deqn{E[T \mid T > t_c] = t_c + \frac{1}{S(t_c)} \int_{t_c}^{T_{\max}} S(t)\,dt}
#' Event subjects are left unchanged.
#'
#' @param time Numeric vector of observed survival/censoring times.
#' @param status Integer or logical vector; 1 (or TRUE) = event, 0 (or FALSE) = censored.
#'
#' @return Numeric vector the same length as \code{time} with censored values
#'   replaced by their KM conditional mean estimates.
#'
#' @examples
#' \dontrun{
#' time   <- c(5, 10, 3, 8, 12)
#' status <- c(1,  0, 1, 0,  1)
#' kmImputeTimes(time, status)
#' }
#'
#' @importFrom survival Surv survfit
#' @export
kmImputeTimes <- function(time, status) {
  km.fit   <- survival::survfit(survival::Surv(time, status) ~ 1)
  km.times <- km.fit$time
  km.surv  <- km.fit$surv
  imputed  <- time

  for (i in which(status == 0)) {
    tc         <- time[i]
    idx.le     <- which(km.times <= tc)
    S.tc       <- if (length(idx.le) == 0) 1.0 else km.surv[max(idx.le)]
    if (S.tc <= 0) next

    idx.beyond <- which(km.times > tc)
    if (length(idx.beyond) == 0) next

    t.knots    <- c(tc, km.times[idx.beyond])
    s.knots    <- c(S.tc, km.surv[idx.beyond])
    area       <- sum(diff(t.knots) * s.knots[-length(s.knots)])
    imputed[i] <- tc + area / S.tc
  }
  imputed
}
