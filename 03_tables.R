###############################################################################
# 03_tables.R
# STAGE 3 -- TABLES
#
# Table 1   Descriptive statistics (H1 sample; primary onset sample by onset)
# Table 1b  Adult-onset counts: offending, illicit drug use, or both
# Table 2a  Correlations, H1 sample
# Table 2b  Correlations, primary onset sample (incl. abuse and ADHD)
# Table 3   H1 models (Wave I delinquency)
# Table 4   Primary adult-onset models (Waves III-V): A, B, C, D
# Table 5   Offending onset vs. illicit drug use onset (Model C)
# Table 6   Robustness: secondary design (all adult onset) and the earlier
#           Wave I-only design
# Table S1  Sample flow;  Table S2  Weight plan
#
# Each table is saved as CSV and all together in output/tables/ALL_TABLES.html
###############################################################################

cat("\n================ STAGE 3: TABLES ================\n")
dat  <- readRDS(file.path(work_dir, "analysis_data.rds"))
mr   <- readRDS(file.path(work_dir, "model_results.rds"))
flow <- read.csv(file.path(log_dir, "sample_flow.csv"))
ct   <- mr$coef_table

fmt   <- function(x, k = 2) formatC(x, format = "f", digits = k)
stars <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
w_mean <- function(x, w) { ok <- !is.na(x); sum(x[ok] * w[ok]) / sum(w[ok]) }

# Analytic samples (exactly the cases each main model used)
h1  <- dat[dat$in_h1, ]
h2p <- dat[dat$in_h2p, ]
h2s <- dat[dat$in_h2s, ]

# ---------------------------------------------------------------------------
# TABLE 1: DESCRIPTIVES
# ---------------------------------------------------------------------------
desc_vars <- c(
  delinq_w1         = "Any Wave I delinquency (%)",
  delinq_variety_w1 = "Wave I delinquency variety (0-13)",
  strain_raw        = "Strain: CES-D (0-57)",
  peer_raw          = "Peer deviance (0-9)",
  bond_raw          = "Bond deficit (mean z; higher = weaker)",
  phys_abuse        = "Physical abuse before 18 (%)",
  sex_abuse         = "Sexual abuse before 18 (%)",
  adhd_count        = "Childhood ADHD symptoms (0-17)",
  adhd_probable     = "Probable childhood ADHD (%)",
  w1_binge          = "Wave I binge drinking, past year (%)",
  age_w1            = "Age at Wave I",
  male              = "Male (%)",
  hisp              = "Hispanic (%)",
  black             = "Non-Hispanic Black (%)",
  other_race        = "Non-Hispanic other race (%)",
  par_college       = "Resident parent has college degree (%)",
  par_ed_unknown    = "Parent education unknown (%)"
)
binary <- c("delinq_w1", "male", "hisp", "black", "other_race", "par_college",
            "par_ed_unknown", "phys_abuse", "sex_abuse", "adhd_probable", "w1_binge")
describe <- function(d, wt) {
  sapply(names(desc_vars), function(v) {
    x <- d[[v]]; w <- d[[wt]]
    if (all(is.na(x))) return("--")
    if (v %in% binary) sprintf("%s / %s", fmt(100 * w_mean(x, w), 1), fmt(100 * mean(x, na.rm = TRUE), 1))
    else sprintf("%s / %s (%s)", fmt(w_mean(x, w)), fmt(mean(x, na.rm = TRUE)), fmt(sd(x, na.rm = TRUE)))
  })
}
t1 <- data.frame(
  Variable = unname(desc_vars),
  H1_sample           = describe(h1, "GSWGT1"),
  Primary_all         = describe(h2p, "GSW1345"),
  Primary_no_onset    = describe(h2p[h2p$onset_all == 0, ], "GSW1345"),
  Primary_onset       = describe(h2p[h2p$onset_all == 1, ], "GSW1345"),
  Secondary_all       = describe(h2s, "GSW1345"),
  Secondary_no_onset  = describe(h2s[h2s$onset_late == 0, ], "GSW1345"),
  Secondary_onset     = describe(h2s[h2s$onset_late == 1, ], "GSW1345"),
  row.names = NULL)
t1[1:2, 3:8] <- "--"
t1 <- rbind(t1, data.frame(Variable = "Unweighted N", H1_sample = nrow(h1),
  Primary_all = nrow(h2p), Primary_no_onset = sum(h2p$onset_all == 0),
  Primary_onset = sum(h2p$onset_all == 1), Secondary_all = nrow(h2s),
  Secondary_no_onset = sum(h2s$onset_late == 0), Secondary_onset = sum(h2s$onset_late == 1)))
t1_note <- paste("Weighted mean / unweighted mean (unweighted SD) for continuous variables;",
  "weighted % / unweighted % for binary variables. Weights: GSWGT1 (H1), GSW1345 (onset).",
  "Abuse and ADHD are measured only for Wave III/IV respondents, so the H1 column describes",
  "those H1 cases who have them. Onset samples are Wave I abstainers from delinquency and illicit",
  "drugs (Wave I delinquency = 0 by design).")

# ---------------------------------------------------------------------------
# TABLE 1b: ONSET COUNTS BY ROUTE
# ---------------------------------------------------------------------------
count_rows <- function(d, wt, any, off, drug, hard, sample_label) {
  tot <- sum(d[[wt]])
  one <- function(label, x) { ok <- !is.na(x)
    data.frame(Sample = sample_label, Indicator = label, n = sum(x[ok] == 1),
               pct_weighted = fmt(100 * sum(d[[wt]][ok] * x[ok]) / sum(d[[wt]][ok]), 1)) }
  r <- d[[sub("onset", "route", any)]]
  rbind(one("ANY ADULT ONSET", d[[any]]),
        one("  Self-reported offending (any)", d[[off]]),
        one("  Illicit drug use (any)", d[[drug]]),
        if (!is.null(hard)) one("    ... excluding marijuana", d[[hard]]),
        one("  Offending only", as.numeric(r == "offending_only")),
        one("  Illicit drug use only", as.numeric(r == "substance_only")),
        one("  Both", as.numeric(r == "both")))
}
t1b <- rbind(count_rows(h2p, "GSW1345", "onset_all", "off_all", "sub_all", "hard_all",
                        "Primary: all adult onset (W3-W5)"),
             count_rows(h2s, "GSW1345", "onset_late", "off_late", "sub_late", NULL,
                        "Secondary: late onset only (W4-W5)"))
t1b_note <- paste("Illicit drug use = marijuana, cocaine, crystal meth, heroin, other illegal drugs,",
  "or nonmedical prescription drugs (Wave V); alcohol is not included. Arrest is not used.",
  "'Offending only', 'drug use only', and 'both' split the onset cases with complete data on both components.")

# ---------------------------------------------------------------------------
# TABLE 2: CORRELATIONS
# ---------------------------------------------------------------------------
cor_table <- function(d, vars, labs, wt) {
  m <- as.matrix(d[, vars]); ok <- complete.cases(m); m <- m[ok, ]; w <- d[[wt]][ok]
  rw <- cov2cor(cov.wt(m, wt = w)$cov); ru <- cor(m); k <- length(vars)
  out <- matrix("", k, k, dimnames = list(labs, as.character(seq_len(k))))
  for (i in 1:k) for (j in 1:k) {
    if (i == j) out[i, j] <- "1"
    if (i > j)  out[i, j] <- fmt(rw[i, j])
    if (i < j)  out[i, j] <- paste0(fmt(ru[i, j]), stars(cor.test(m[, i], m[, j])$p.value))
  }
  data.frame(Variable = paste0(seq_len(k), ". ", labs), out, row.names = NULL, check.names = FALSE)
}
t2a <- cor_table(h1, c("delinq_w1", "strain_raw", "peer_raw", "bond_raw", "age_w1", "male", "par_college"),
                 c("Wave I delinquency", "Strain", "Peer deviance", "Bond deficit", "Age", "Male", "Parent college"),
                 "GSWGT1")
t2b <- cor_table(h2p, c("onset_all", "strain_raw", "peer_raw", "bond_raw", "phys_abuse", "sex_abuse",
                        "adhd_count", "w1_binge", "age_w1", "male", "par_college"),
                 c("Any adult onset", "Strain", "Peer deviance", "Bond deficit", "Physical abuse",
                   "Sexual abuse", "ADHD symptoms", "W1 binge drinking", "Age", "Male", "Parent college"),
                 "GSW1345")
t2_note <- paste("Below the diagonal: weighted Pearson correlations (primary). Above: unweighted;",
  "* p<.05 ** p<.01 *** p<.001 (unweighted tests). Binary-variable correlations are point-biserial/phi.")

# ---------------------------------------------------------------------------
# MODEL TABLES
# ---------------------------------------------------------------------------
term_labels <- c(strain = "Strain (CES-D, z)", peer = "Peer deviance (z)",
                 bond = "Bond deficit (z)", "strain:peer" = "Strain x Peer deviance",
                 phys_abuse = "Physical abuse before 18", sex_abuse = "Sexual abuse before 18",
                 adhd = "Childhood ADHD symptoms (z)", w1_binge = "Wave I binge drinking",
                 age_w1 = "Age at Wave I", male = "Male", hisp = "Hispanic",
                 black = "Non-Hispanic Black", other_race = "Non-Hispanic other",
                 par_college = "Parent college degree", par_ed_unknown = "Parent education unknown",
                 "(Intercept)" = "Intercept")
model_col <- function(spec, model, weighting = "weighted", family = "logit") {
  x <- ct[ct$spec == spec & ct$model == model & ct$weighting == weighting & ct$family == family, ]
  cell <- if (family == "logit") sprintf("%s [%s, %s]%s", fmt(x$OR), fmt(x$OR_lo), fmt(x$OR_hi), stars(x$p))
          else sprintf("%s (%s)%s", fmt(x$b, 3), fmt(x$se, 3), stars(x$p))
  list(cells = setNames(cell, x$term), n = x$n[1], ev = x$n_events[1])
}
build_model_table <- function(cols) {
  terms <- names(term_labels)
  out <- data.frame(Term = unname(term_labels), row.names = NULL)
  for (nm in names(cols)) out[[nm]] <- ifelse(terms %in% names(cols[[nm]]$cells), cols[[nm]]$cells[terms], "")
  out <- out[apply(out[, -1, drop = FALSE], 1, function(r) any(r != "")), ]
  rbind(out, c("N", sapply(cols, `[[`, "n")), c("Onset / outcome events", sapply(cols, `[[`, "ev")))
}
model_note <- paste("Logistic regression odds ratios [95% CI]; LPM = linear probability model,",
  "b (SE) = change in probability. Weighted models use the weight plan (Table S2) with CLUSTER2",
  "cluster-robust SEs; unweighted models keep the same cluster correction.",
  "Strain, peer deviance, bond deficit, and ADHD symptoms are z-scores.",
  "* p<.05 ** p<.01 *** p<.001.")

t3 <- build_model_table(list(
  "Model A weighted OR"   = model_col("H1", "A"),
  "Model B weighted OR"   = model_col("H1", "B"),
  "Model B unweighted OR" = model_col("H1", "B", "unweighted"),
  "Model B weighted LPM"  = model_col("H1", "B", family = "LPM")))

t4 <- build_model_table(list(
  "A weighted OR"   = model_col("H2P", "A"),
  "B weighted OR"   = model_col("H2P", "B"),
  "C weighted OR"   = model_col("H2P", "C"),
  "C unweighted OR" = model_col("H2P", "C", "unweighted"),
  "C weighted LPM"  = model_col("H2P", "C", family = "LPM"),
  "D (+ W1 binge drinking) weighted OR" = model_col("H2P", "D")))

t5 <- build_model_table(list(
  "Primary (W3-W5): offending"                   = model_col("H2P_OFF", "C"),
  "Primary (W3-W5): illicit drug use"            = model_col("H2P_SUB", "C"),
  "Primary: drug use excl. marijuana"    = model_col("H2P_HARD", "C"),
  "Late onset: offending"                 = model_col("H2S_OFF", "C"),
  "Late onset: illicit drug use"          = model_col("H2S_SUB", "C")))
t5_note <- paste(model_note, "Model C, weighted. Each column is a separate outcome estimated on the full",
  "risk set (onset of that kind = 1, otherwise 0). Columns with few events (bottom row) are imprecise.")

t6 <- build_model_table(list(
  "Late onset (W3 screen) B weighted" = model_col("H2S", "B"),
  "Late onset C weighted"                   = model_col("H2S", "C"),
  "Late onset C weighted LPM"               = model_col("H2S", "C", family = "LPM"),
  "Earlier W1-only design B weighted"      = model_col("H2X", "B"),
  "Earlier W1-only design B unweighted"    = model_col("H2X", "B", "unweighted")))

# ---------------------------------------------------------------------------
# SAVE
# ---------------------------------------------------------------------------
tables <- list(
  "Table 1. Descriptive statistics" = list(t1, t1_note),
  "Table 1b. Adult-onset counts: offending and illicit drug use" = list(t1b, t1b_note),
  "Table 2a. Correlations, H1 sample" = list(t2a, t2_note),
  "Table 2b. Correlations, primary onset sample" = list(t2b, t2_note),
  "Table 3. Wave I delinquency models (H1)" = list(t3, model_note),
  "Table 4. Adult-onset models (H2 primary: Waves III-V)" = list(t4, model_note),
  "Table 5. Offending onset vs. illicit drug use onset" = list(t5, t5_note),
  "Table 6. Robustness: late-onset design (W3 screen) and the Wave I-only design" = list(t6, model_note),
  "Table S1. Sample flow (N after each filter)" = list(flow, "N after every filtering step, from the Wave I public-use file."),
  "Table S2. Weight plan" = list(weight_plan, "Decided before any model was run (see 00_setup.R).")
)
file_names <- c("table1_descriptives", "table1b_onset_counts", "table2a_correlations_H1",
                "table2b_correlations_onset", "table3_models_H1", "table4_models_onset_primary",
                "table5_models_offending_vs_drug", "table6_robustness", "tableS1_sample_flow",
                "tableS2_weight_plan")
for (i in seq_along(tables)) write.csv(tables[[i]][[1]], file.path(tab_dir, paste0(file_names[i], ".csv")), row.names = FALSE)

html_table <- function(df) {
  head <- paste0("<tr>", paste0("<th>", names(df), "</th>", collapse = ""), "</tr>")
  rows <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table>", head, paste(rows, collapse = ""), "</table>")
}
html <- c("<html><head><meta charset='utf-8'><title>Add Health adult-onset tables</title>",
  "<style>body{font-family:Georgia,serif;margin:30px;max-width:1400px}",
  "table{border-collapse:collapse;margin-bottom:6px;font-size:13px}",
  "th,td{border:1px solid #bbb;padding:4px 8px;text-align:left}",
  "th{background:#eee}p.note{font-size:12px;color:#444;margin-bottom:34px}</style></head><body>",
  paste0("<h1>Add Health adult-onset analysis: tables</h1><p>Generated ", Sys.time(), "</p>"))
for (nm in names(tables)) html <- c(html, paste0("<h2>", nm, "</h2>"), html_table(tables[[nm]][[1]]),
                                    paste0("<p class='note'>", tables[[nm]][[2]], "</p>"))
writeLines(c(html, "</body></html>"), file.path(tab_dir, "ALL_TABLES.html"))
cat("Saved", length(tables), "tables to output/tables/ (CSV + ALL_TABLES.html)\n")
print(t1b, row.names = FALSE)
