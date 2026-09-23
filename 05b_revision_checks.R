###############################################################################
# 05b_revision_checks.R
# STAGE 5b -- THREE CHECKS A REVIEWER WILL ASK FOR (detection paper)
#
#   CHECK 1  Does perceived peer substance use still predict first arrest once
#            the respondent's OWN adolescent substance use is controlled?
#            (The alternative reading is projection: people who use report
#            that their friends use.)
#   CHECK 2  Are the composition contrasts statistically significant?
#            Adult-onset arrestees vs the never-arrested, share with juvenile
#            offending, under the broad and strict screens. Reported as a
#            weighted prevalence ratio with a 95% CI, plus a Rao-Scott test.
#   CHECK 3  Do the composition shares hold WITHOUT the Wave IV-V first
#            arrests whose ages were imputed?
#
# Runs after 05_arrest_onset.R (reads data_work/arrest_onset_results.rds).
# Writes output/arrest_onset/revision_checks.csv and prints a summary.
###############################################################################

cat("\n================ STAGE 5b: REVISION CHECKS ================\n")
ao <- readRDS(file.path(work_dir, "arrest_onset_results.rds"))
d <- ao$d; py1 <- ao$py1; py2 <- ao$py2; w3_block <- ao$w3_block
fmt <- function(x, k = 2) formatC(x, format = "f", digits = k)
stars <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))

# Same model setup as Stage 5
w1_block  <- c("strain", "peer", "bond", "parent_control_z", controls)
adv_block <- adversity
rhs_age   <- "age_c + I(age_c^2)"
fit_py <- function(py, rhs) {
  des <- svydesign(ids = ~CLUSTER2, weights = ~GSW1345, data = py)
  svyglm(as.formula(paste("event ~", rhs)), design = des, family = quasibinomial())
}
or_row <- function(fit, term) {
  b <- unname(coef(fit)[term]); se <- unname(sqrt(diag(vcov(fit)))[term])
  p <- 2 * pt(-abs(b / se), fit$df.residual)
  c(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = p)
}
cell <- function(r) sprintf("%s [%s, %s]%s", fmt(r["OR"]), fmt(r["lo"]), fmt(r["hi"]), stars(r["p"]))

out <- list()

# ---------------------------------------------------------------------------
# CHECK 1: OWN ADOLESCENT SUBSTANCE USE
# ---------------------------------------------------------------------------
# w1_su = z-score of the 0-3 count of Wave I substance use types: past-year
# binge drinking, ever marijuana, ever cocaine/inhalants/other illegal drugs.
# (The arrest risk set is NOT screened on Wave I drug use, so the full range
# is present.) A second version enters the three types separately.
cat("\n-- CHECK 1: peer substance use, net of own adolescent substance use --\n")
for (w in list(list("Window 1", py1, c(w1_block, adv_block)),
               list("Window 2", py2, c(w1_block, adv_block, w3_block)))) {
  py <- w[[2]][complete.cases(w[[2]][, c("w1_su", "w1_binge", "w1_marij", "w1_hard")]), ]
  base <- paste(rhs_age, "+", paste(w[[3]], collapse = " + "))
  f0 <- fit_py(py, base)
  f1 <- fit_py(py, paste(base, "+ w1_su"))
  f2 <- fit_py(py, paste(base, "+ w1_binge + w1_marij + w1_hard"))
  rows <- list(
    c(w[[1]], "Peer substance use, same sample, no own-use control", cell(or_row(f0, "peer"))),
    c(w[[1]], "Peer substance use, + own use (count, z)",            cell(or_row(f1, "peer"))),
    c(w[[1]], "  Own adolescent substance use (count, z)",           cell(or_row(f1, "w1_su"))),
    c(w[[1]], "Peer substance use, + own use (three types)",         cell(or_row(f2, "peer"))),
    c(w[[1]], "  Own past-year binge drinking",                      cell(or_row(f2, "w1_binge"))),
    c(w[[1]], "  Own ever marijuana",                                cell(or_row(f2, "w1_marij"))),
    c(w[[1]], "  Own ever other illegal drugs",                      cell(or_row(f2, "w1_hard"))))
  ev <- sum(py$event)
  for (r in rows) { out[[length(out) + 1]] <- c("Check 1", r, ev)
                    cat(sprintf("  %-9s %-55s %s\n", r[1], r[2], r[3])) }
  cat(sprintf("  %-9s people = %d, first arrests = %d\n", w[[1]],
              length(unique(py$AID)), sum(py$event)))
}

# ---------------------------------------------------------------------------
# CHECK 2 AND 3: COMPOSITION CONTRASTS, WITH AND WITHOUT IMPUTED AGES
# ---------------------------------------------------------------------------
# Prevalence ratio = weighted share with juvenile offending among adult-onset
# arrestees / weighted share among the never-arrested. Estimated with a
# log-link quasi-Poisson survey model (robust SEs), which gives the ratio and
# its CI directly. The Rao-Scott chi-square is reported alongside.
rs <- d[d$in_win1, ]
rs$adult_onset <- as.numeric(!is.na(rs$age_first))
contrast <- function(dd, screen, label) {
  dd <- dd[!is.na(dd[[screen]]), ]
  des <- svydesign(ids = ~CLUSTER2, weights = ~GSW1345, data = dd)
  f <- svyglm(as.formula(paste(screen, "~ adult_onset")), design = des, family = quasipoisson())
  b <- coef(f)["adult_onset"]; se <- sqrt(vcov(f)["adult_onset", "adult_onset"])
  chi <- svychisq(as.formula(paste("~", screen, "+ adult_onset")), design = des, statistic = "F")
  w <- dd$GSW1345; g <- dd$adult_onset; y <- dd[[screen]]
  p_on <- 100 * sum(w[g == 1] * y[g == 1]) / sum(w[g == 1])
  p_no <- 100 * sum(w[g == 0] * y[g == 0]) / sum(w[g == 0])
  r <- c(label, fmt(p_on, 1), fmt(p_no, 1),
         sprintf("%s [%s, %s]", fmt(exp(b)), fmt(exp(b - 1.96 * se)), fmt(exp(b + 1.96 * se))),
         ifelse(chi$p.value < .001, "<.001", fmt(chi$p.value, 3)),
         sum(g == 1), sum(g == 0))
  cat(sprintf("  %-46s onset %5s%%  never %5s%%  PR %-20s p %s  (n %s / %s)\n",
              r[1], r[2], r[3], r[4], r[5], r[6], r[7]))
  r
}
cat("\n-- CHECKS 2-3: juvenile offending, adult-onset arrestees vs never-arrested --\n")
no_imp <- rs[!(rs$adult_onset == 1 & rs$age_first_imputed %in% TRUE), ]
comp <- list(
  contrast(rs,     "sr_juv_delinq", "Broad screen, all adult-onset arrests"),
  contrast(rs,     "sr_juv_strict", "Strict screen, all adult-onset arrests"),
  contrast(no_imp, "sr_juv_delinq", "Broad screen, reported ages only"),
  contrast(no_imp, "sr_juv_strict", "Strict screen, reported ages only"))
comp_df <- do.call(rbind, comp)
colnames(comp_df) <- c("Contrast", "Adult_onset_pct", "Never_arrested_pct",
                       "Prevalence_ratio_95CI", "Rao_Scott_p", "n_onset", "n_never")

# ---------------------------------------------------------------------------
# SAVE
# ---------------------------------------------------------------------------
c1 <- do.call(rbind, out)
colnames(c1) <- c("Check", "Window", "Row", "OR_95CI", "First_arrests")
write.csv(c1, file.path(out_dir, "arrest_onset", "revision_check1_own_substance_use.csv"), row.names = FALSE)
write.csv(comp_df, file.path(out_dir, "arrest_onset", "revision_checks2_3_composition.csv"), row.names = FALSE)
saveRDS(list(check1 = c1, composition = comp_df), file.path(work_dir, "revision_checks.rds"))
cat("\nSaved output/arrest_onset/revision_check1_own_substance_use.csv\n")
cat("Saved output/arrest_onset/revision_checks2_3_composition.csv\n")
