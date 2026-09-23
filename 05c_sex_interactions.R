###############################################################################
# 05c_sex_interactions.R
# STAGE 5c -- DOES ANY OF THIS DIFFER BY SEX?
#
# Ragan and Jacobsen (2026, Criminology) report that the consequences of
# juvenile arrest for peer relationships are stronger for females. Reviewers
# from that literature will ask whether the ANTECEDENTS of first arrest are
# likewise gendered. Sex is only a control in the main models, so this script
# adds the test:
#
#   (1) Interactions of male with each focal predictor, one model per window.
#   (2) Sex-stratified models, so the associations can be read directly for
#       women and for men.
#   (3) The same for the Wave III adult roles in Window 2.
#
# Runs after 05_arrest_onset.R. Writes output/arrest_onset/sex_interactions.csv
# and prints a summary.
###############################################################################

cat("\n================ STAGE 5c: SEX DIFFERENCES ================\n")
ao <- readRDS(file.path(work_dir, "arrest_onset_results.rds"))
py1 <- ao$py1; py2 <- ao$py2; w3_block <- ao$w3_block; lab <- ao$lab
fmt <- function(x, k = 2) formatC(x, format = "f", digits = k)
stars <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))

w1_block  <- c("strain", "peer", "bond", "parent_control_z", controls)
adv_block <- adversity
rhs_age   <- "age_c + I(age_c^2)"

fit_py <- function(py, rhs) {
  des <- svydesign(ids = ~CLUSTER2, weights = ~GSW1345, data = py)
  svyglm(as.formula(paste("event ~", rhs)), design = des, family = quasibinomial())
}
row_of <- function(fit, term, model, window, n_people, n_events) {
  b <- coef(fit); se <- sqrt(diag(vcov(fit)))
  alt <- if (grepl(":", term)) paste(rev(strsplit(term, ":")[[1]]), collapse = ":") else term
  k <- intersect(c(term, alt), names(b))
  if (!length(k)) return(NULL)
  k <- k[1]; p <- 2 * pt(-abs(b[k] / se[k]), fit$df.residual)
  data.frame(Window = window, Model = model, Term = term,
             OR = sprintf("%s [%s, %s]%s", fmt(exp(unname(b[k]))),
                          fmt(exp(unname(b[k] - 1.96 * se[k]))),
                          fmt(exp(unname(b[k] + 1.96 * se[k]))), stars(unname(p))),
             People = n_people, Events = n_events, row.names = NULL)
}

focal <- c("peer", "strain", "bond", "parent_control_z", "phys_abuse", "sex_abuse", "adhd")
out <- list()

for (w in list(list("Window 1", py1, c(w1_block, adv_block), focal),
               list("Window 2", py2, c(w1_block, adv_block, w3_block), c(focal, w3_block)))) {
  win <- w[[1]]; py <- w[[2]]; blk <- w[[3]]; terms <- w[[4]]
  base <- paste(rhs_age, "+", paste(blk, collapse = " + "))

  # (1) male x predictor, entered one at a time so that each test is clean
  for (t in terms) {
    f <- fit_py(py, paste(base, "+ male:", t))
    out[[length(out) + 1]] <- row_of(f, paste0("male:", t), "male x predictor", win,
                                     length(unique(py$AID)), sum(py$event))
  }
  # (2) stratified by sex
  for (g in c(0, 1)) {
    pg <- py[py$male == g, ]
    gl <- ifelse(g == 1, "men only", "women only")
    blk_g <- setdiff(blk, "male")
    f <- fit_py(pg, paste(rhs_age, "+", paste(blk_g, collapse = " + ")))
    for (t in terms) out[[length(out) + 1]] <-
      row_of(f, t, gl, win, length(unique(pg$AID)), sum(pg$event))
  }
}
res <- do.call(rbind, out)
res$Term <- ifelse(res$Term %in% names(lab), unlist(lab[res$Term]),
                   sub("^male:", "Male x ", res$Term))
write.csv(res, file.path(out_dir, "arrest_onset", "sex_interactions.csv"), row.names = FALSE)

cat("\n-- Interactions of male with each predictor (one at a time) --\n")
int <- res[res$Model == "male x predictor", ]
for (i in seq_len(nrow(int)))
  cat(sprintf("  %-9s %-34s %s\n", int$Window[i], int$Term[i], int$OR[i]))
cat("\n-- Stratified models (read the focal rows against each other) --\n")
st <- res[res$Model != "male x predictor", ]
for (i in seq_len(nrow(st)))
  cat(sprintf("  %-9s %-11s %-34s %s  (events %s)\n", st$Window[i], st$Model[i],
              st$Term[i], st$OR[i], st$Events[i]))
cat("\nSaved output/arrest_onset/sex_interactions.csv\n")
cat("Note: the stratified models are the more readable of the two, but the\n",
    "interaction terms are the formal test. Report the interaction p values\n",
    "rather than comparing significance across the stratified columns.\n")
