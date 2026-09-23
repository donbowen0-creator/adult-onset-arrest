###############################################################################
# 06_manuscript_tables.R
# STAGE 6 -- MANUSCRIPT TABLES AND FIGURES (detection paper, JCJ)
#
# Builds the journal-ready tables and figure files from the Stage 5 results:
#   Table 1  Sample descriptives by adult-onset arrest status
#   Table 2  Composition of adult-onset arrest: juvenile offending under two
#            screens, and arrest charge types
#   Table 3  Discrete-time hazard models, window 1 and window 2
#   Table 4  Robustness of the peer substance use effect
#   Table 5  Pooled windows: predictor x window interactions
#   Figure 1 Age-specific hazard of first arrest (reported ages only)
#   Figure 2 Peer substance use vs. co-offending and non-substance deviance
#
# Outputs go to output/manuscript/ :
#   table1_descriptives.csv ... table5_pooled.csv   (read by make_tables_docx.js)
#   Figure_1.pdf / Figure_1.png, Figure_2.pdf / Figure_2.png
#   (JCJ wants each figure as its own file; PDF is the vector version.)
###############################################################################

cat("\n================ STAGE 6: MANUSCRIPT TABLES ================\n")
ao  <- readRDS(file.path(work_dir, "arrest_onset_results.rds"))
d <- ao$d; coefs <- ao$coefs; lab <- ao$lab; py1 <- ao$py1; py2 <- ao$py2
man_dir <- file.path(out_dir, "manuscript"); dir.create(man_dir, showWarnings = FALSE)

fmt   <- function(x, k = 2) formatC(x, format = "f", digits = k)
stars <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
wm <- function(x, w) { ok <- !is.na(x); sum(x[ok] * w[ok]) / sum(w[ok]) }

# ---------------------------------------------------------------------------
# TABLE 1: DESCRIPTIVES (window 1 risk set; roles on the window 2 risk set)
# ---------------------------------------------------------------------------
rs <- d[d$in_win1, ]
rs$adult_onset <- as.numeric(!is.na(rs$age_first))
desc_rows <- list(
  list("Perceived peer substance use (0-9)", "peer_raw", FALSE),
  list("Depressive symptoms, CES-D (0-57)", "strain_raw", FALSE),
  list("Bond deficit (mean z)", "bond_raw", FALSE),
  list("Parental control of decisions (0-7)", "parent_control", FALSE),
  list("Physical abuse before 18 (%)", "phys_abuse", TRUE),
  list("Sexual abuse before 18 (%)", "sex_abuse", TRUE),
  list("Childhood ADHD symptoms (0-17)", "adhd_count", FALSE),
  list("Any juvenile delinquency, broad screen (%)", "sr_juv_delinq", TRUE),
  list("Serious juvenile delinquency, strict screen (%)", "sr_juv_strict", TRUE),
  list("Age at Wave I", "age_w1", FALSE),
  list("Male (%)", "male", TRUE),
  list("Hispanic (%)", "hisp", TRUE),
  list("Non-Hispanic Black (%)", "black", TRUE),
  list("Non-Hispanic other race (%)", "other_race", TRUE),
  list("Resident parent has college degree (%)", "par_college", TRUE))
cell <- function(x, w, pct) {
  if (pct) fmt(100 * wm(x, w), 1)
  else sprintf("%s (%s)", fmt(wm(x, w)), fmt(sd(x, na.rm = TRUE)))
}
t1 <- do.call(rbind, lapply(desc_rows, function(r) {
  x <- rs[[r[[2]]]]; w <- rs$GSW1345; g <- rs$adult_onset
  data.frame(Variable = r[[1]], Total = cell(x, w, r[[3]]),
             Never_arrested = cell(x[g == 0], w[g == 0], r[[3]]),
             Adult_onset = cell(x[g == 1], w[g == 1], r[[3]]))
}))
# Wave III roles, window 2 risk set
r2 <- d[d$in_win2, ]; r2$adult_onset <- as.numeric(!is.na(r2$fail_w2))
role_rows <- list(list("Out of school (%)", "school_exit_w3"), list("Left parental home (%)", "left_home_w3"),
                  list("Working 10+ hours/week (%)", "working_w3"), list("Married or cohabiting (%)", "partner_w3"),
                  list("Children in household (%)", "children_w3"), list("Material hardship (%)", "hardship_w3"))
t1b <- do.call(rbind, lapply(role_rows, function(r) {
  x <- r2[[r[[2]]]]; w <- r2$GSW1345; g <- r2$adult_onset
  data.frame(Variable = paste0("  ", r[[1]]), Total = cell(x, w, TRUE),
             Never_arrested = cell(x[g == 0], w[g == 0], TRUE),
             Adult_onset = cell(x[g == 1], w[g == 1], TRUE))
}))
t1 <- rbind(t1, data.frame(Variable = "Wave III roles (window 2 risk set)", Total = "", Never_arrested = "", Adult_onset = ""),
            t1b,
            data.frame(Variable = "Unweighted n (window 1 risk set)", Total = nrow(rs),
                       Never_arrested = sum(rs$adult_onset == 0), Adult_onset = sum(rs$adult_onset == 1)),
            data.frame(Variable = "Unweighted n (window 2 risk set)", Total = nrow(r2),
                       Never_arrested = sum(r2$adult_onset == 0), Adult_onset = sum(r2$adult_onset == 1)))

# ---------------------------------------------------------------------------
# TABLE 2: COMPOSITION -- juvenile offending, and arrest charge types
# ---------------------------------------------------------------------------
paths <- locate_icpsr()
w4c <- as.data.frame(read_dta(paths$w4_main, col_select = c("AID", paste0("H4CJ7", LETTERS[1:11]))))
ch4 <- function(v) { x <- clean_num(w4c[[v]], c(6, 8, 9)); x[x == 7] <- 0; x }
w4c$sub4  <- as.numeric(rowSums(sapply(paste0("H4CJ7", LETTERS[1:4]),  ch4) == 1, na.rm = TRUE) > 0)
w4c$dui4  <- as.numeric(rowSums(sapply(paste0("H4CJ7", LETTERS[1:2]),  ch4) == 1, na.rm = TRUE) > 0)
w4c$drug4 <- as.numeric(rowSums(sapply(paste0("H4CJ7", LETTERS[3:4]),  ch4) == 1, na.rm = TRUE) > 0)
w4c$oth4  <- as.numeric(rowSums(sapply(paste0("H4CJ7", LETTERS[5:11]), ch4) == 1, na.rm = TRUE) > 0)
w4c$any4  <- as.numeric(w4c$sub4 == 1 | w4c$oth4 == 1)
w5o <- paste0("H5CJ4", c("E", "F", "H", "I", "J", "K", "L"))
w5c <- as.data.frame(read_dta(paths$w5_main, col_select = all_of(c("AID", paste0("H5CJ4", LETTERS[1:4]), w5o))))
ch5 <- function(v) { x <- clean_num(w5c[[v]], 97); x[is.na(x)] <- 0; x }
w5c$sub5  <- as.numeric(rowSums(sapply(paste0("H5CJ4", LETTERS[1:4]), ch5) == 1) > 0)
w5c$dui5  <- as.numeric(rowSums(sapply(paste0("H5CJ4", LETTERS[1:2]), ch5) == 1) > 0)
w5c$drug5 <- as.numeric(rowSums(sapply(paste0("H5CJ4", LETTERS[3:4]), ch5) == 1) > 0)
w5c$oth5  <- as.numeric(rowSums(sapply(w5o, ch5) == 1) > 0)
w5c$any5  <- as.numeric(w5c$sub5 == 1 | w5c$oth5 == 1)
ao_c <- merge(rs[rs$adult_onset == 1, c("AID", "GSW1345", "age_first_imputed")],
              w4c[, c("AID", "sub4", "dui4", "drug4", "oth4", "any4")], by = "AID", all.x = TRUE)
ao_c <- merge(ao_c, w5c[, c("AID", "sub5", "dui5", "drug5", "oth5", "any5")], by = "AID", all.x = TRUE)
# Arrests dated by Wave IV use Wave IV charges; Wave IV-V arrests use Wave V charges.
pick <- function(a, b) ifelse(ao_c$age_first_imputed, ao_c[[b]], ao_c[[a]])
ao_c$sub <- pick("sub4", "sub5"); ao_c$dui <- pick("dui4", "dui5"); ao_c$drug <- pick("drug4", "drug5")
ao_c$oth <- pick("oth4", "oth5"); ao_c$anych <- pick("any4", "any5")
ao_c$subonly <- as.numeric(ao_c$sub == 1 & ao_c$oth == 0)
known <- ao_c[ao_c$anych %in% 1, ]
cat(sprintf("Adult-onset arrestees with a charge type reported: %d of %d\n", nrow(known), nrow(ao_c)))
pc <- function(v) fmt(100 * wm(known[[v]], known$GSW1345), 1)

noarr <- rs[rs$adult_onset == 0, ]; onset <- rs[rs$adult_onset == 1, ]
t2 <- data.frame(
  Measure = c("A. Juvenile self-reported offending (Waves I-II)",
              "  Any delinquency, broad screen (%)", "  Serious delinquency, strict screen (%)",
              "B. Charges for adult-onset arrests (among those reporting a charge)",
              "  Any substance-related charge (%)", "    DUI/DWI or other alcohol (%)",
              "    Marijuana or other drug (%)", "  Substance-related charges only (%)",
              "  Any non-substance charge (%)", "Unweighted n"),
  Never_arrested = c("", fmt(100 * wm(noarr$sr_juv_delinq, noarr$GSW1345), 1),
                     fmt(100 * wm(noarr$sr_juv_strict, noarr$GSW1345), 1),
                     "", "--", "--", "--", "--", "--", nrow(noarr)),
  Adult_onset = c("", fmt(100 * wm(onset$sr_juv_delinq, onset$GSW1345), 1),
                  fmt(100 * wm(onset$sr_juv_strict, onset$GSW1345), 1),
                  "", pc("sub"), pc("dui"), pc("drug"), pc("subonly"), pc("oth"),
                  sprintf("%d (charges: %d)", nrow(onset), nrow(known))))
# Prevalence ratios and Rao-Scott tests from Stage 5b
rc <- readRDS(file.path(work_dir, "revision_checks.rds"))
cmp <- rc$composition
pr_row <- function(label, i) data.frame(Measure = label, Never_arrested = "(ref.)",
  Adult_onset = sprintf("%s, p %s", cmp[i, "Prevalence_ratio_95CI"],
                        ifelse(cmp[i, "Rao_Scott_p"] == "<.001", "< .001", paste("=", cmp[i, "Rao_Scott_p"]))))
t2 <- rbind(t2[1:3, ],
            pr_row("  Prevalence ratio, broad screen [95% CI]", 1),
            pr_row("  Prevalence ratio, strict screen [95% CI]", 2),
            t2[4:nrow(t2), ])
# Unweighted counts too, since the table rows above are weighted
cat(sprintf("Unweighted: broad %d/%d, strict %d/%d among adult-onset arrestees\n",
            sum(onset$sr_juv_delinq %in% 1), nrow(onset), sum(onset$sr_juv_strict %in% 1), nrow(onset)))

# ---------------------------------------------------------------------------
# TABLES 3-5: MODELS
# ---------------------------------------------------------------------------
ci <- function(model, term) {
  x <- coefs[coefs$model == model & coefs$term == term, ]
  if (!nrow(x)) {
    alt <- paste(rev(strsplit(term, ":")[[1]]), collapse = ":")
    x <- coefs[coefs$model == model & coefs$term == alt, ]
  }
  if (!nrow(x)) return("")
  sprintf("%s [%s, %s]%s", fmt(x$OR), fmt(x$OR_lo), fmt(x$OR_hi), stars(x$p))
}
meta <- function(model, what) coefs[[what]][coefs$model == model][1]
t3_terms <- c("peer", "strain", "bond", "parent_control_z", "phys_abuse", "sex_abuse", "adhd",
              "school_exit_w3", "left_home_w3", "working_w3", "partner_w3", "children_w3", "hardship_w3",
              "age_c", "I(age_c^2)", "age_w1", "male", "hisp", "black", "other_race",
              "par_college", "par_ed_unknown")
m1 <- "Window 1: Wave I + adversity"; m2 <- "Window 2: + Wave III roles"
t3 <- data.frame(Predictor = sapply(t3_terms, function(t) lab[[t]]),
                 Window_1 = sapply(t3_terms, function(t) ci(m1, t)),
                 Window_2 = sapply(t3_terms, function(t) ci(m2, t)), row.names = NULL)
t3$Predictor[t3$Predictor == "Peer deviance (z)"] <- "Perceived peer substance use (z)"
t3$Predictor[t3$Predictor == "Strain (CES-D, z)"] <- "Depressive symptoms, CES-D (z)"
t3 <- rbind(t3, data.frame(Predictor = c("Persons", "Person-years", "First arrests"),
  Window_1 = c(meta(m1, "people"), meta(m1, "person_years"), meta(m1, "events")),
  Window_2 = c(meta(m2, "people"), meta(m2, "person_years"), meta(m2, "events"))))

rob <- list(
  c("Baseline", "Window 1: Wave I + adversity", "Window 2: + Wave III roles", "peer"),
  c("Controlling own reported offending", "Window 1: + own Wave I delinquency variety",
    "Window 2: + own Wave III offending", "peer"),
  c("No self-reported juvenile delinquency (broad screen)", "Window 1: no self-reported juvenile delinquency",
    "Window 2: no self-reported juvenile delinquency", "peer"),
  c("No serious juvenile delinquency (strict screen)", "Window 1: no SERIOUS juvenile delinquency",
    "Window 2: no SERIOUS juvenile delinquency", "peer"),
  c("Network subsample, all peer measures entered: perceived peer substance use",
    "Window 1 split: all four", "Window 2 split: all four", "peer"),
  c("  Friends' own substance use", "Window 1 split: all four", "Window 2 split: all four", "fr_su_z"),
  c("  Friends' own non-substance risk behavior", "Window 1 split: all four", "Window 2 split: all four", "fr_risk_z"),
  c("  Own group offending", "Window 1 split: all four", "Window 2 split: all four", "group_fight"))
t4 <- do.call(rbind, lapply(rob, function(r) data.frame(
  Specification = r[1],
  Window_1 = sprintf("%s (%s)", ci(r[2], r[4]), meta(r[2], "events")),
  Window_2 = sprintf("%s (%s)", ci(r[3], r[4]), meta(r[3], "events")))))
c1 <- as.data.frame(rc$check1, stringsAsFactors = FALSE)
own <- function(w) { k <- c1$Window == w & c1$Row == "Peer substance use, + own use (count, z)"
  sprintf("%s (%s)", c1$OR_95CI[k], c1$First_arrests[k]) }
t4 <- rbind(t4[1:2, ],
            data.frame(Specification = "Controlling own adolescent substance use",
                       Window_1 = own("Window 1"), Window_2 = own("Window 2")),
            t4[3:nrow(t4), ])
rep_m <- "Reported ages only: 18 to Wave IV"
t4 <- rbind(t4, data.frame(Specification = "Reported ages only, age 18 to Wave IV (pooled)",
  Window_1 = sprintf("%s (%s)", ci(rep_m, "peer"), meta(rep_m, "events")), Window_2 = "--"))

pm <- "Pooled windows: predictor x window"
t5_terms <- c("late:peer", "late:strain", "late:phys_abuse", "late:sex_abuse", "late:adhd")
t5 <- data.frame(Interaction = sapply(t5_terms, function(t) lab[[t]]),
                 OR = sapply(t5_terms, function(t) ci(pm, t)), row.names = NULL)
t5$Interaction <- sub("Peer deviance", "Perceived peer substance use", t5$Interaction)
t5$Interaction <- sub("^Strain", "Depressive symptoms", t5$Interaction)
t5 <- rbind(t5, data.frame(Interaction = c("Persons", "First arrests"),
                           OR = c(meta(pm, "people"), meta(pm, "events"))))

# ---------------------------------------------------------------------------
# TABLE 6: SEX DIFFERENCES (built from Stage 5c)
# ---------------------------------------------------------------------------
sx_file <- file.path(out_dir, "arrest_onset", "sex_interactions.csv")
if (file.exists(sx_file)) {
  sx <- read.csv(sx_file, stringsAsFactors = FALSE)
  pick <- function(win, mod, term) {
    v <- sx$OR[sx$Window == win & sx$Model == mod & sx$Term == term]
    if (length(v)) v[1] else ""
  }
  terms6 <- unique(sx$Term[sx$Model == "male x predictor" & sx$Window == "Window 2"])
  raw6 <- sub("^Male x ", "", terms6)
  pretty6 <- sapply(raw6, function(r) if (r %in% names(lab)) lab[[r]] else r)
  t6 <- data.frame(
    Predictor = pretty6,
    W1_interaction = sapply(terms6, function(t) pick("Window 1", "male x predictor", t)),
    W2_interaction = sapply(terms6, function(t) pick("Window 2", "male x predictor", t)),
    W2_women = sapply(pretty6, function(t) pick("Window 2", "women only", t)),
    W2_men = sapply(pretty6, function(t) pick("Window 2", "men only", t)),
    row.names = NULL)
  write.csv(t6, file.path(man_dir, "table6_sex_differences.csv"), row.names = FALSE)
}

write.csv(t1, file.path(man_dir, "table1_descriptives.csv"), row.names = FALSE)
write.csv(t2, file.path(man_dir, "table2_composition.csv"), row.names = FALSE)
write.csv(t3, file.path(man_dir, "table3_hazard_models.csv"), row.names = FALSE)
write.csv(t4, file.path(man_dir, "table4_peer_robustness.csv"), row.names = FALSE)
write.csv(t5, file.path(man_dir, "table5_pooled_windows.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# FIGURES (each its own file, as JCJ requires)
# ---------------------------------------------------------------------------
theme_set(theme_classic(base_size = 11, base_family = "serif"))
drep <- d[d$in_win1, ]
drep$fail_rep <- ifelse(!drep$age_first_imputed, drep$age_first, NA)
# Rebuild the reported-age person-years (same rules as Stage 5)
pyr <- do.call(rbind, lapply(seq_len(nrow(drep)), function(i) {
  a1 <- floor(drep$age_w4[i]); f <- drep$fail_rep[i]
  if (!is.finite(a1) || a1 < 18) return(NULL)
  if (!is.na(f) && f >= 18 && f <= a1) a1 <- floor(f)
  ages <- 18:a1
  data.frame(age = ages, GSW1345 = drep$GSW1345[i],
             event = as.numeric(!is.na(f) & f >= 18 & ages == floor(f)))
}))
hz <- aggregate(cbind(e = event * GSW1345, w = GSW1345) ~ age, data = pyr, sum)
hz$py <- as.vector(table(pyr$age)[as.character(hz$age)])
hz <- hz[hz$py >= 500, ]; hz$hazard <- 100 * hz$e / hz$w
f1 <- ggplot(hz, aes(age, hazard)) + geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_y_continuous(limits = c(0, NA)) + scale_x_continuous(breaks = seq(18, 30, 2)) +
  labs(x = "Age", y = "First arrests per 100 at risk")
ggsave(file.path(man_dir, "Figure_1.pdf"), f1, width = 6, height = 3.8)
ggsave(file.path(man_dir, "Figure_1.png"), f1, width = 6, height = 3.8, dpi = 600, bg = "white")

sp <- c(peer = "Perceived peer substance use", fr_su_z = "Friends' own substance use",
        fr_risk_z = "Friends' own non-substance risk", group_fight = "Own group offending")
f2d <- do.call(rbind, lapply(c("Window 1", "Window 2"), function(w) {
  x <- coefs[coefs$model == paste(w, "split: all four") & coefs$term %in% names(sp), ]
  data.frame(window = paste0(w, ifelse(w == "Window 1", ": age 18 to Wave III", ": Wave III onward")),
             term = factor(sp[x$term], levels = rev(sp)), OR = x$OR, lo = x$OR_lo, hi = x$OR_hi)
}))
f2 <- ggplot(f2d, aes(OR, term)) + geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.35) + facet_wrap(~ window) +
  scale_x_log10(breaks = c(0.7, 0.8, 1, 1.25, 1.5, 2)) + labs(x = "Odds ratio for first arrest (log scale), 95% CI", y = NULL) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(man_dir, "Figure_2.pdf"), f2, width = 7, height = 3.2)
ggsave(file.path(man_dir, "Figure_2.png"), f2, width = 7, height = 3.2, dpi = 600, bg = "white")

cat("\n-- Table 2B: charges for adult-onset arrests (weighted %) --\n")
print(t2[5:9, c("Measure", "Adult_onset")], row.names = FALSE)
cat("\nSaved manuscript tables (CSV) and Figure_1 / Figure_2 (PDF + PNG) to output/manuscript/\n")
