###############################################################################
# 02_models.R
# STAGE 2 -- MODELS
#
# Nested logistic models:
#   Model A = strain + peer deviance + bond deficit + controls
#   Model B = A + strain x peer deviance
#   Model C = B + physical abuse + sexual abuse + childhood ADHD symptoms
#   Model D = C + own Wave I binge drinking            (robustness)
#
# Outcomes / samples (see 00_setup.R, sections 3g-3i):
#   H1        Wave I delinquency                          A, B        GSWGT1
#   H2P       primary: any adult onset (W3 window)         A, B, C, D  GSW1345
#   H2P_OFF   primary: self-reported offending onset      C           GSW1345
#   H2P_SUB   primary: illicit drug use onset             C           GSW1345
#   H2P_HARD  primary: drug use onset, marijuana excluded C           GSW1345
#   H2S       secondary: late onset only (W3 screen)       A, B, C     GSW1345
#   H2S_OFF   secondary: offending onset                  C           GSW1345
#   H2S_SUB   secondary: illicit drug use onset           C           GSW1345
#   H2X       sensitivity: W1-only screen                 A, B        GSW145
#
# Every model is run WEIGHTED (primary) and UNWEIGHTED; both use
# cluster-robust SEs (CLUSTER2), so the only difference is the weight.
# The top model of each main outcome is also re-run as a weighted linear
# probability model (LPM) to check the interaction is not a logit artifact.
#
# Weighted designs are built on EVERYONE who has the weight, and the analytic
# cases are selected with subset() -- the correct subpopulation analysis.
#
# Produces: data_work/model_results.rds, output/tables/all_coefficients.csv
###############################################################################

cat("\n================ STAGE 2: MODELS ================\n")
dat <- readRDS(file.path(work_dir, "analysis_data.rds"))

rhs <- list()
rhs$A <- paste(c("strain", "peer", "bond", controls), collapse = " + ")
rhs$B <- paste(rhs$A, "+ strain:peer")
rhs$C <- paste(rhs$B, "+", paste(adversity, collapse = " + "))
rhs$D <- paste(rhs$C, "+ w1_binge")

sp <- function(label, y, flag, wt, models, lpm = NULL)
  list(label = label, y = y, flag = flag, wt = wt, models = models, lpm = lpm)
specs <- list(
  H1       = sp("H1: Wave I delinquency", "delinq_w1", "in_h1", "GSWGT1", c("A", "B"), "B"),
  H2P      = sp("H2 primary: any adult onset (W3-W5)", "onset_all", "in_h2p", "GSW1345", c("A", "B", "C", "D"), "C"),
  H2P_OFF  = sp("H2 primary: offending onset", "off_all", "in_h2p", "GSW1345", "C"),
  H2P_SUB  = sp("H2 primary: illicit drug use onset", "sub_all", "in_h2p", "GSW1345", "C"),
  H2P_HARD = sp("H2 primary: drug use onset excl. marijuana", "hard_all", "in_h2p", "GSW1345", "C"),
  H2S      = sp("H2 secondary: late onset only (W4-W5)", "onset_late", "in_h2s", "GSW1345", c("A", "B", "C"), "C"),
  H2S_OFF  = sp("H2 secondary: offending onset", "off_late", "in_h2s", "GSW1345", "C"),
  H2S_SUB  = sp("H2 secondary: illicit drug use onset", "sub_late", "in_h2s", "GSW1345", "C"),
  H2X      = sp("H2 sensitivity: W1-only screen", "onset_w45", "in_h2x", "GSW145", c("A", "B"))
)

tidy_fit <- function(fit, spec_id, model, weighting, n, n_yes, family = "logit") {
  est <- coef(fit); se <- sqrt(diag(vcov(fit)))
  p <- 2 * pt(-abs(est / se), df = fit$df.residual)
  out <- data.frame(spec = spec_id, model = model, weighting = weighting,
                    family = family, term = names(est), b = est, se = se,
                    p = p, row.names = NULL)
  if (family == "logit") {
    out$OR <- exp(out$b); out$OR_lo <- exp(out$b - 1.96 * out$se); out$OR_hi <- exp(out$b + 1.96 * out$se)
  } else out$OR <- out$OR_lo <- out$OR_hi <- NA
  out$n <- n; out$n_events <- n_yes
  out
}

fits <- list(); coefs <- list()
for (id in names(specs)) {
  sp <- specs[[id]]
  cat("\n----", sp$label, "----\n")
  d <- dat[!is.na(dat[[sp$wt]]), ]
  # Route-specific outcomes are NA for other routes -> those cases leave the sample
  d$analytic <- d[[sp$flag]] & !is.na(d[[sp$y]])
  needed <- unique(c(sp$y, all.vars(as.formula(paste("~", rhs[[tail(sp$models, 1)]])))))
  d$analytic <- d$analytic & complete.cases(d[, needed])
  # Out-of-sample rows keep zero domain weight; fill their blanks so R does
  # not silently drop them from the cluster structure.
  for (v in needed) d[[v]][!d$analytic & is.na(d[[v]])] <- 0
  n_a <- sum(d$analytic); n_y <- sum(d[[sp$y]][d$analytic])
  cat(sprintf("  N = %d, events = %d\n", n_a, n_y))

  des_w <- subset(svydesign(ids = ~CLUSTER2, weights = as.formula(paste0("~", sp$wt)), data = d), analytic)
  des_u <- svydesign(ids = ~CLUSTER2, weights = ~1, data = d[d$analytic, ])

  for (m in sp$models) {
    f <- as.formula(paste(sp$y, "~", rhs[[m]]))
    fw <- svyglm(f, design = des_w, family = quasibinomial())
    fu <- svyglm(f, design = des_u, family = quasibinomial())
    fits[[paste(id, m, "weighted", sep = "_")]]   <- fw
    fits[[paste(id, m, "unweighted", sep = "_")]] <- fu
    coefs[[length(coefs) + 1]] <- tidy_fit(fw, id, m, "weighted", n_a, n_y)
    coefs[[length(coefs) + 1]] <- tidy_fit(fu, id, m, "unweighted", n_a, n_y)
  }
  if (!is.null(sp$lpm)) {
    fl <- svyglm(as.formula(paste(sp$y, "~", rhs[[sp$lpm]])), design = des_w, family = gaussian())
    fits[[paste(id, sp$lpm, "weighted_LPM", sep = "_")]] <- fl
    coefs[[length(coefs) + 1]] <- tidy_fit(fl, id, sp$lpm, "weighted", n_a, n_y, family = "LPM")
  }
  key <- do.call(rbind, coefs)
  key <- key[key$spec == id & key$weighting == "weighted" &
             key$term %in% c("strain", "peer", "bond", "strain:peer", adversity, "w1_binge"),
             c("model", "family", "term", "b", "se", "p")]
  key$OR <- round(exp(key$b), 2); key$OR[key$family == "LPM"] <- NA
  key[, c("b", "se")] <- round(key[, c("b", "se")], 3); key$p <- signif(key$p, 3)
  print(key, row.names = FALSE)
  fits[[paste(id, "design_w", sep = "_")]] <- des_w
}

coef_table <- do.call(rbind, coefs)
write.csv(coef_table, file.path(tab_dir, "all_coefficients.csv"), row.names = FALSE)
saveRDS(list(fits = fits, specs = specs, coef_table = coef_table, rhs = rhs),
        file.path(work_dir, "model_results.rds"))
cat("\nSaved data_work/model_results.rds and output/tables/all_coefficients.csv\n")
