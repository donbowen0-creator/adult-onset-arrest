###############################################################################
# 05_arrest_onset.R
# STAGE 5 -- ADULT-ONSET ARREST: TWO STACKED RISK SETS (discrete-time survival)
#
# Outcome: age at FIRST arrest (Wave IV: H4CJ2/H4CJ3/H4CJ4), updated with
# Wave V for anyone with no arrest by Wave IV who reports one at Wave V
# (H5CJ3; Wave V records no age, so the failure is placed at the midpoint of
# the Wave IV and Wave V interview ages and flagged).
#   at risk from age 18 | failure = first arrest at 18+ |
#   first arrest before 18 = excluded (already failed) |
#   never arrested by the last interview = censored at that age.
#
# WINDOW 1 -- "expiry": person-years from age 18 to the Wave III interview
#   age. Predictors are Wave I adolescent strain, bonds, and scaffold
#   (parental decision-making, living with parents) plus childhood adversity.
#   Wave III role measures are NOT used here: they are recorded at the END of
#   this window, so for a failure at 20 with Wave III at 22 they would
#   postdate the event.
#
# WINDOW 2 -- "late": person-years from the Wave III interview age to the
#   last interview age, for people still arrest-free at Wave III. Wave III
#   roles (school exit, left home, work, partner, children in household,
#   material hardship) are fully lagged predictors of a later first arrest.
#   Role QUALITY is tested with work x hardship and partner x hardship
#   interactions: Sampson-Laub predicts protective main effects; adult GST
#   predicts the interaction carries it.
#
# SENSITIVITY (Beckley/Moffitt): official adult onset split by whether the
#   person self-reported any delinquency at Wave I or Wave II. Wave II is
#   used ONLY for this split.
#
# Output: output/arrest_onset/  (CSV tables, ARREST_ONSET.html, 3 figures)
###############################################################################

cat("\n================ STAGE 5: ADULT-ONSET ARREST ================\n")
paths <- locate_icpsr()
dat <- readRDS(file.path(work_dir, "analysis_data.rds"))
ao_dir <- file.path(out_dir, "arrest_onset"); dir.create(ao_dir, showWarnings = FALSE)
flow_log <- flow_log[0, ]

# ---------------------------------------------------------------------------
# 1. AGES AND AGE AT FIRST ARREST
# ---------------------------------------------------------------------------
w1a <- as.data.frame(read_dta(paths$w1_main, col_select = c("AID", "H1GI1M", "H1GI1Y",
        paste0("H1WP", 1:7))))
# ICPSR stores some years two-digit (94) and some four-digit (2016); put
# every year on the same footing before doing any date arithmetic.
yr4 <- function(y) ifelse(is.na(y), NA, ifelse(y < 30, y + 2000, ifelse(y < 100, y + 1900, y)))
bm <- clean_num(w1a$H1GI1M, 96); by <- yr4(clean_num(w1a$H1GI1Y, 96))
w1a$birth <- by * 12 + bm                       # birth month, absolute scale
# Adolescent scaffold: number of the 7 decisions the parents still control
# (low autonomy = tight scaffold), 0-7, needs >= 5 answered
auton <- sapply(paste0("H1WP", 1:7), function(v) clean_num(w1a[[v]], c(6, 7, 8, 9)))
n_a <- rowSums(!is.na(auton))
w1a$parent_control <- ifelse(n_a >= 5, rowSums(1 - auton, na.rm = TRUE) * 7 / n_a, NA)

w4a <- as.data.frame(read_dta(paths$w4_main, col_select = c("AID", "IMONTH4", "IYEAR4",
        "H4CJ1", "H4CJ2", "H4CJ3", "H4CJ4")))
cj1 <- clean_num(w4a$H4CJ1, c(6, 7, 8)); cj2 <- clean_num(w4a$H4CJ2, c(6, 7, 8))
cj3 <- clean_num(w4a$H4CJ3, c(96, 97, 98)); cj4 <- clean_num(w4a$H4CJ4, c(96, 97, 98))
w4a$ever_arr_w4 <- cj1
w4a$age_first_arrest <- ifelse(cj1 %in% 0, NA, ifelse(cj2 %in% 1, cj3, cj4))
w4a$int4 <- yr4(clean_num(w4a$IYEAR4, NA)) * 12 + clean_num(w4a$IMONTH4, NA)

w5a <- as.data.frame(read_dta(paths$w5_main, col_select = c("AID", "IMONTH5", "IYEAR5", "H5CJ3")))
w5a$int5 <- yr4(clean_num(w5a$IYEAR5, NA)) * 12 + clean_num(w5a$IMONTH5, NA)
w5a$ever_arr_w5 <- clean_num(w5a$H5CJ3, numeric(0))

d <- merge(dat, w1a[, c("AID", "birth", "parent_control")], by = "AID", all.x = TRUE)
d <- merge(d, w4a[, c("AID", "int4", "ever_arr_w4", "age_first_arrest")], by = "AID", all.x = TRUE)
d <- merge(d, w5a[, c("AID", "int5", "ever_arr_w5")], by = "AID", all.x = TRUE)
d$age_w4 <- (d$int4 - d$birth) / 12
d$age_w5 <- (d$int5 - d$birth) / 12

# Wave V update: first-ever arrest between Waves IV and V (age not recorded)
d$arrest_w4_w5 <- ifelse(d$ever_arr_w4 %in% 0 & d$ever_arr_w5 %in% 1, 1, 0)
d$age_first_imputed <- d$arrest_w4_w5 == 1
d$age_first <- ifelse(d$arrest_w4_w5 == 1, (d$age_w4 + d$age_w5) / 2, d$age_first_arrest)
d$censor_age <- ifelse(!is.na(d$age_w5), d$age_w5, d$age_w4)

cat("\nAge at first arrest, Wave I + III + IV + V respondents with GSW1345:\n")
base <- !is.na(d$GSW1345) & d$in_w3 == 1 & d$in_w4 == 1 & d$in_w5 == 1
print(table(cut(d$age_first[base], c(-Inf, 17.99, 24.99, Inf),
                labels = c("first arrest < 18", "18-24", "25+")),
            ifelse(d$age_first_imputed[base], "age imputed (W4-W5)", "age reported"),
            useNA = "ifany"))

# ---------------------------------------------------------------------------
# 2. WAVE III ROLES (lagged predictors for window 2)
# ---------------------------------------------------------------------------
hard_items <- paste0("H3EC", 18:24)
w3r <- as.data.frame(read_dta(paths$w3_main, col_select = all_of(c("AID", "H3ED23",
        "H3HR2", "H3LM7", "H3MR1", "H3MR13_A", "H3DA18", "H3DA19", "H3CJ1", "H3CJ3",
        hard_items))))
# Arrested by the Wave III interview (used to break ties when the reported
# age at first arrest equals the Wave III age: the arrest could fall either
# side of the interview, and Wave III says which).
stopped3 <- clean_num(w3r$H3CJ1, c(96, 98, 99, 996, 998, 999))
arr3 <- clean_num(w3r$H3CJ3, c(6, 8, 9))
arr3[arr3 == 7 & stopped3 %in% 0] <- 0; arr3[arr3 == 7] <- NA
w3r$arrest_by_w3 <- arr3
w3r$in_school_w3 <- clean_num(w3r$H3ED23, c(8, 9))
w3r$school_exit_w3 <- 1 - w3r$in_school_w3
res <- clean_num(w3r$H3HR2, c(6, 9, 99))
w3r$left_home_w3 <- ifelse(res %in% 1, 0, ifelse(res %in% 2:5, 1, NA))
wk <- clean_num(w3r$H3LM7, c(6, 8, 9)); wk[wk == 7] <- NA
w3r$working_w3 <- wk
mar <- clean_num(w3r$H3MR1, c(6, 8, 9))
coh <- clean_num(w3r$H3MR13_A, c(6, 8, 9)); coh[coh == 7] <- 0
w3r$partner_w3 <- ifelse(mar >= 1 | coh %in% 1, 1, ifelse(!is.na(mar), 0, NA))
kids <- sapply(c("H3DA18", "H3DA19"), function(v) clean_num(w3r[[v]], c(96, 98, 99)))
w3r$children_w3 <- as.numeric(rowSums(kids > 0, na.rm = TRUE) > 0)
hard <- sapply(hard_items, function(v) { x <- clean_num(w3r[[v]], c(6, 8, 9)); x[x == 7] <- 0; x })
w3r$hardship_w3 <- any_act(as.data.frame(hard), colnames(hard))
w3r <- w3r[, c("AID", "arrest_by_w3", "school_exit_w3", "left_home_w3", "working_w3",
               "partner_w3", "children_w3", "hardship_w3")]
d <- merge(d, w3r, by = "AID", all.x = TRUE)

# ---------------------------------------------------------------------------
# 3. WAVE II DELINQUENCY (Beckley/Moffitt sensitivity split only)
# ---------------------------------------------------------------------------
w2_items <- c("H2DS1", "H2DS2", "H2DS4", "H2DS6", "H2DS7", "H2DS8", "H2DS9",
              "H2DS10", "H2DS11", "H2DS13", "H2FV6", "H2FV16", "H2FV22")
w2_path <- file.path(raw_dir, "ICPSR_21600/DS0005/21600-0005-Data.dta")
if (!file.exists(w2_path)) {
  zips <- c(if (file.exists(icpsr_source) && !dir.exists(icpsr_source)) icpsr_source,
            list.files(raw_dir, pattern = "^ICPSR_21600.*\\.zip$", full.names = TRUE))
  if (length(zips)) unzip(zips[1], files = "ICPSR_21600/DS0005/21600-0005-Data.dta", exdir = raw_dir)
}
w2 <- as.data.frame(read_dta(w2_path, col_select = all_of(c("AID", w2_items))))
for (v in w2_items) w2[[v]] <- clean_num(w2[[v]], c(6, 8, 9))
w2$H2FV22[w2$H2FV22 == 7] <- 0
w2$delinq_w2 <- any_act(w2, w2_items)
# STRICT screen: only acts no one would call adolescent nuisance -- serious
# theft, burglary, car theft, weapon use, drug sales, injuring someone.
w2_strict <- c("H2DS6", "H2DS7", "H2DS8", "H2DS9", "H2DS10", "H2FV6", "H2FV22")
w2$delinq_w2_strict <- any_act(w2, w2_strict)
d <- merge(d, w2[, c("AID", "delinq_w2", "delinq_w2_strict")], by = "AID", all.x = TRUE)
w1_strict <- c("H1DS6", "H1DS8", "H1DS9", "H1DS10", "H1DS11", "H1DS12", "H1FV7")
w1s <- as.data.frame(read_dta(paths$w1_main, col_select = all_of(c("AID", w1_strict))))
for (v in w1_strict) w1s[[v]] <- clean_num(w1s[[v]], c(6, 8, 9))
w1s$delinq_w1_strict <- any_act(w1s, w1_strict)
d <- merge(d, w1s[, c("AID", "delinq_w1_strict")], by = "AID", all.x = TRUE)
# Any self-reported juvenile delinquency at Wave I or Wave II
d$sr_juv_delinq <- ifelse(d$delinq_w1 %in% 1 | d$delinq_w2 %in% 1, 1,
                   ifelse(d$delinq_w1 %in% 0, 0, NA))
# Same thing under the strict screen (serious acts only)
d$sr_juv_strict <- ifelse(d$delinq_w1_strict %in% 1 | d$delinq_w2_strict %in% 1, 1,
                   ifelse(d$delinq_w1_strict %in% 0, 0, NA))

# ---------------------------------------------------------------------------
# 3b. SPLITTING PEER DEVIANCE: SUBSTANCE USE vs GROUP OFFENDING
# ---------------------------------------------------------------------------
# The main peer-deviance scale (H1TO9/29/33) is entirely FRIENDS' SUBSTANCE
# USE as the respondent reports it. Add Health public-use has no direct
# co-offending measure, so the split uses the closest available pieces:
#   fr_su    -- friends' OWN reports of smoking, drinking, getting drunk
#               (in-school survey items S59A-C, averaged over nominated
#               friends; Wave I network file). Substance, but not filtered
#               through the respondent's perception.
#   fr_risk  -- friends' OWN reports of non-substance risk behaviour
#               (S59D-G: racing, dangerous dares, lying to parents, skipping
#               school). Deviance, not offending.
#   group_fight -- the respondent's OWN group offending at Wave I (H1DS14:
#               fight where a group of friends was against another group).
#               The one co-offending act in the Wave I instrument.
# Item wording for S59A-G is from the In-School Questionnaire, which is not
# in the public zip; confirm it against the Add Health documentation.
net <- as.data.frame(read_dta(paths$w1_net, col_select = c("AID",
         paste0("AXS59", LETTERS[1:7]))))
net$fr_su   <- rowMeans(net[, paste0("AXS59", LETTERS[1:3])], na.rm = TRUE)
net$fr_risk <- rowMeans(net[, paste0("AXS59", LETTERS[4:7])], na.rm = TRUE)
for (v in c("fr_su", "fr_risk")) net[[v]][is.nan(net[[v]])] <- NA
net$fr_su_z <- zscore(net$fr_su); net$fr_risk_z <- zscore(net$fr_risk)
w1g <- as.data.frame(read_dta(paths$w1_main, col_select = c("AID", "H1DS14")))
w1g$group_fight <- as.numeric(clean_num(w1g$H1DS14, c(6, 8, 9)) > 0)
d <- merge(d, net[, c("AID", "fr_su_z", "fr_risk_z")], by = "AID", all.x = TRUE)
d <- merge(d, w1g[, c("AID", "group_fight")], by = "AID", all.x = TRUE)

# ---------------------------------------------------------------------------
# 4. RISK SETS AND SAMPLE FLOW
# ---------------------------------------------------------------------------
w1_block  <- c("strain", "peer", "bond", "parent_control_z", controls)
adv_block <- adversity
w3_block  <- c("school_exit_w3", "left_home_w3", "working_w3", "partner_w3",
               "children_w3", "hardship_w3")
d$parent_control_z <- zscore(d$parent_control)

cat("\n---- Sample flow: WINDOW 1 (expiry: age 18 to Wave III) ----\n")
s <- rep(TRUE, nrow(d))
record_step("W1 expiry", "Wave I public-use in-home respondents", sum(s))
s <- s & d$in_w3 == 1 & d$in_w4 == 1 & d$in_w5 == 1
record_step("W1 expiry", "interviewed at Waves III, IV and V", sum(s))
s <- s & !is.na(d$GSW1345);  record_step("W1 expiry", "has longitudinal weight (GSW1345)", sum(s))
s <- s & !(d$age_first %in% 0:17.99) & !(d$age_first < 18 & !is.na(d$age_first))
record_step("W1 expiry", "no first arrest before age 18", sum(s))
s <- s & !is.na(d$age_w3) & d$age_w3 > 18
record_step("W1 expiry", "older than 18 at Wave III (has exposure)", sum(s))
s <- s & complete.cases(d[, c(w1_block, adv_block)])
record_step("W1 expiry", "complete Wave I + adversity predictors = RISK SET", sum(s))
d$in_win1 <- s

# Split the failures between the windows. A first arrest at an age below the
# Wave III age belongs to window 1; above it, to window 2. When the reported
# age EQUALS the Wave III age, the Wave III arrest report decides.
before_w3 <- !is.na(d$age_first) &
  (d$age_first < d$age_w3 | (d$age_first == d$age_w3 & d$arrest_by_w3 %in% 1))
d$fail_w1 <- ifelse(before_w3, d$age_first, NA)
d$fail_w2 <- ifelse(!is.na(d$age_first) & !before_w3, d$age_first, NA)

cat("\n---- Sample flow: WINDOW 2 (late: Wave III to last interview) ----\n")
s <- d$in_win1 & is.na(d$fail_w1)
record_step("W2 late", "window 1 risk set, still arrest-free at Wave III", sum(s))
s <- s & d$censor_age > d$age_w3
record_step("W2 late", "observed after Wave III", sum(s))
s <- s & complete.cases(d[, w3_block])
record_step("W2 late", "complete Wave III role measures = RISK SET", sum(s))
d$in_win2 <- s

# ---------------------------------------------------------------------------
# 5. PERSON-YEAR FILES
# ---------------------------------------------------------------------------
# One row per person per year of exposure; the outcome is 1 in the year the
# first arrest happens and 0 in every earlier year. Age is kept as a
# predictor, so this is a discrete-time hazard model.
make_py <- function(dd, start, stop, failcol = "age_first") {
  out <- lapply(seq_len(nrow(dd)), function(i) {
    a0 <- floor(start[i]); a1 <- floor(min(stop[i], dd$censor_age[i], na.rm = TRUE))
    fail_age <- dd[[failcol]][i]
    failed <- !is.na(fail_age) && fail_age >= a0 && fail_age <= a1
    if (failed) a1 <- floor(fail_age)
    if (!is.finite(a0) || !is.finite(a1) || a1 < a0 || a1 - a0 > 30) return(NULL)
    ages <- a0:a1
    data.frame(AID = dd$AID[i], age = ages,
               event = as.numeric(failed & ages == max(ages)))
  })
  py <- do.call(rbind, out)
  merge(py, dd, by = "AID")
}
py1 <- make_py(d[d$in_win1, ], rep(18, sum(d$in_win1)), d$age_w3[d$in_win1], "fail_w1")
py2 <- make_py(d[d$in_win2, ], d$age_w3[d$in_win2], d$censor_age[d$in_win2], "fail_w2")
py1$age_c <- py1$age - 18
py2$age_c <- py2$age - 22

cat(sprintf("\nWindow 1: %d people, %d person-years, %d first arrests (%.1f%% of people)\n",
            sum(d$in_win1), nrow(py1), sum(py1$event), 100 * sum(py1$event) / sum(d$in_win1)))
cat(sprintf("Window 2: %d people, %d person-years, %d first arrests (%.1f%% of people)\n",
            sum(d$in_win2), nrow(py2), sum(py2$event), 100 * sum(py2$event) / sum(d$in_win2)))

# ---------------------------------------------------------------------------
# 6. MODELS (discrete-time logistic hazard, weighted, cluster-robust)
# ---------------------------------------------------------------------------
fit_py <- function(py, rhs, label) {
  des <- svydesign(ids = ~CLUSTER2, weights = ~GSW1345, data = py)
  f <- svyglm(as.formula(paste("event ~", rhs)), design = des, family = quasibinomial())
  b <- coef(f); se <- sqrt(diag(vcov(f))); p <- 2 * pt(-abs(b / se), f$df.residual)
  data.frame(model = label, term = names(b), b = b, se = se, p = p,
             OR = exp(b), OR_lo = exp(b - 1.96 * se), OR_hi = exp(b + 1.96 * se),
             people = length(unique(py$AID)), person_years = nrow(py),
             events = sum(py$event), row.names = NULL)
}
rhs_age1 <- "age_c + I(age_c^2)"
rhs_age2 <- "age_c + I(age_c^2)"

res <- list()
# Window 1: adolescent strain/bonds/scaffold + childhood adversity
res$W1_core <- fit_py(py1, paste(rhs_age1, "+", paste(c(w1_block, adv_block), collapse = " + ")),
                      "Window 1: Wave I + adversity")
# (The strain x peer interaction is the companion paper's question and is
#  deliberately not estimated here.)
# Window 2: same block + lagged Wave III roles, then role quality
res$W2_core <- fit_py(py2, paste(rhs_age2, "+", paste(c(w1_block, adv_block), collapse = " + ")),
                      "Window 2: Wave I + adversity")
res$W2_roles <- fit_py(py2, paste(rhs_age2, "+", paste(c(w1_block, adv_block, w3_block), collapse = " + ")),
                       "Window 2: + Wave III roles")
res$W2_qual <- fit_py(py2, paste(rhs_age2, "+", paste(c(w1_block, adv_block, w3_block), collapse = " + "),
                                 "+ working_w3:hardship_w3 + partner_w3:hardship_w3"),
                      "Window 2: + role quality")

# Sensitivity: Beckley/Moffitt split. Failures by people who self-reported
# delinquency at Wave I or II are removed (their person-years are dropped),
# so the outcome is official adult onset with no detected juvenile offending.
py1s <- py1[py1$sr_juv_delinq %in% 0, ]; py2s <- py2[py2$sr_juv_delinq %in% 0, ]
res$W1_clean <- fit_py(py1s, paste(rhs_age1, "+", paste(c(w1_block, adv_block), collapse = " + ")),
                       "Window 1: no self-reported juvenile delinquency")
res$W2_clean <- fit_py(py2s, paste(rhs_age2, "+", paste(c(w1_block, adv_block, w3_block), collapse = " + ")),
                       "Window 2: no self-reported juvenile delinquency")
# Strict-screen version of the same two models
py1t <- py1[py1$sr_juv_strict %in% 0, ]; py2t <- py2[py2$sr_juv_strict %in% 0, ]
res$W1_strict <- fit_py(py1t, paste(rhs_age1, "+", paste(c(w1_block, adv_block), collapse = " + ")),
                        "Window 1: no SERIOUS juvenile delinquency")
res$W2_strict <- fit_py(py2t, paste(rhs_age2, "+", paste(c(w1_block, adv_block, w3_block), collapse = " + ")),
                        "Window 2: no SERIOUS juvenile delinquency")

# REPORTED-AGE FILE: every failure has an age the respondent actually gave.
# Everyone is censored at the Wave IV interview age, which is the last date
# an exact age at first arrest exists for, and the 229 Wave IV-to-Wave V
# arrests with no reported age are dropped rather than placed at a midpoint.
drep <- d[d$in_win1, ]
drep$fail_rep <- ifelse(!drep$age_first_imputed, drep$age_first, NA)
pyr <- make_py(drep, rep(18, nrow(drep)), drep$age_w4, "fail_rep")
pyr <- pyr[pyr$age <= floor(pyr$age_w4), ]
pyr$age_c <- pyr$age - 18
res$REP <- fit_py(pyr, paste("age_c + I(age_c^2) +", paste(c(w1_block, adv_block), collapse = " + ")),
                  "Reported ages only: 18 to Wave IV")
# Does the age gradient differ by childhood self-regulation? A maturation
# account implies a STEEPER decline for those with more ADHD symptoms, i.e.
# a negative adhd x age term. A flat interaction means the age decline is
# common to everyone and does not track self-regulation.
res$REP_AGE <- fit_py(pyr, paste("age_c + I(age_c^2) +", paste(c(w1_block, adv_block), collapse = " + "),
                                 "+ adhd:age_c"), "Reported ages: + ADHD x age")
res$REP_AGE2 <- fit_py(pyr, paste("age_c + I(age_c^2) +", paste(c(w1_block, adv_block), collapse = " + "),
                                  "+ adhd:age_c + peer:age_c"), "Reported ages: + ADHD x age + peer x age")

# DETECTION vs INVOLVEMENT. Two readings of the peer-deviance effect:
#   INVOLVEMENT -- deviant peers mean more offending, so more arrests.
#   DETECTION   -- deviant peers mean more exposure to policed situations and
#                  to group offending, so the same behaviour is caught sooner.
# They separate by conditioning the arrest hazard on the respondent's OWN
# reported offending. If peer deviance survives that conditioning, it is
# adding arrest risk beyond the behaviour the person reports.
res$W1_involve <- fit_py(py1, paste(rhs_age1, "+", paste(c(w1_block, adv_block), collapse = " + "),
                                    "+ delinq_variety_w1"),
                         "Window 1: + own Wave I delinquency variety")
py2$offend_w3_c <- ifelse(is.na(py2$offend_w3), 0, py2$offend_w3)
res$W2_involve <- fit_py(py2, paste(rhs_age2, "+", paste(c(w1_block, adv_block, w3_block), collapse = " + "),
                                    "+ offend_w3_c"),
                         "Window 2: + own Wave III offending")
# The sharpest version: people who reported NO offending at all as juveniles.
# Any peer effect there cannot be their own reported behaviour.
res$W1_nonoff <- fit_py(py1[py1$sr_juv_delinq %in% 0, ],
                        paste(rhs_age1, "+", paste(c(w1_block, adv_block), collapse = " + ")),
                        "Window 1: self-reported non-offenders only")

# SPLIT MODELS, run on the subsample with network data so every column
# shares one N. w1_nopeer / w2_nopeer = the usual blocks minus "peer".
w1_nopeer <- setdiff(c(w1_block, adv_block), "peer")
w2_nopeer <- setdiff(c(w1_block, adv_block, w3_block), "peer")
split_models <- function(py, base, rhs_age, tag) {
  pn <- py[complete.cases(py[, c("fr_su_z", "fr_risk_z", "group_fight")]), ]
  b <- paste(rhs_age, "+", paste(base, collapse = " + "))
  list(
    fit_py(pn, paste(b, "+ peer"), paste(tag, "split: perceived peer substance use")),
    fit_py(pn, paste(b, "+ peer + group_fight"), paste(tag, "split: + own group offending")),
    fit_py(pn, paste(b, "+ fr_su_z + fr_risk_z"), paste(tag, "split: friends' own reports")),
    fit_py(pn, paste(b, "+ peer + fr_su_z + fr_risk_z + group_fight"), paste(tag, "split: all four")))
}
res <- c(res, split_models(py1, w1_nopeer, rhs_age1, "Window 1"),
              split_models(py2, w2_nopeer, rhs_age2, "Window 2"))

# Pooled test: does early adversity act differently in the two windows?
# The windows are disjoint in person-time, so they stack. "late" = window 2.
py1$late <- 0; py2$late <- 1
common <- intersect(names(py1), names(py2))
pyp <- rbind(py1[, common], py2[, common])
pyp$age_c <- pyp$age - 18
res$POOL <- fit_py(pyp, paste("age_c + I(age_c^2) + late +",
  paste(c(w1_block, adv_block), collapse = " + "),
  "+ late:sex_abuse + late:phys_abuse + late:adhd + late:peer + late:strain"),
  "Pooled windows: predictor x window")
coefs <- do.call(rbind, res)

# ---------------------------------------------------------------------------
# 7. TABLES
# ---------------------------------------------------------------------------
fmt <- function(x, k = 2) formatC(x, format = "f", digits = k)
stars <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
lab <- c(late = "Window 2 (vs window 1)",
         fr_su_z = "Friends' own substance use (z)", fr_risk_z = "Friends' own non-substance risk (z)",
         group_fight = "Own group offending (group fight)",
         delinq_variety_w1 = "Own Wave I delinquency variety (0-13)",
         offend_w3_c = "Own Wave III offending (past year)", "age_c:adhd" = "ADHD symptoms x Age",
         "adhd:age_c" = "ADHD symptoms x Age", "age_c:peer" = "Peer deviance x Age",
         "peer:age_c" = "Peer deviance x Age",
         "late:sex_abuse" = "Sexual abuse x Window 2", "late:phys_abuse" = "Physical abuse x Window 2",
         "late:adhd" = "ADHD symptoms x Window 2", "late:peer" = "Peer deviance x Window 2",
         "late:strain" = "Strain x Window 2",
         strain = "Strain (CES-D, z)", peer = "Peer deviance (z)", bond = "Bond deficit (z)",
         parent_control_z = "Parental control of decisions (z)", 
         phys_abuse = "Physical abuse before 18", sex_abuse = "Sexual abuse before 18",
         adhd = "Childhood ADHD symptoms (z)", school_exit_w3 = "Out of school (W3)",
         left_home_w3 = "Left parental home (W3)", working_w3 = "Working 10+ hrs (W3)",
         partner_w3 = "Married or cohabiting (W3)", children_w3 = "Children in household (W3)",
         hardship_w3 = "Material hardship (W3)",
         "working_w3:hardship_w3" = "Working x Hardship", "partner_w3:hardship_w3" = "Partner x Hardship",
         age_c = "Age (centered)", "I(age_c^2)" = "Age squared", age_w1 = "Age at Wave I",
         male = "Male", hisp = "Hispanic", black = "Non-Hispanic Black",
         other_race = "Non-Hispanic other", par_college = "Parent college degree",
         par_ed_unknown = "Parent education unknown", "(Intercept)" = "Intercept")
mk_table <- function(models) {
  terms <- names(lab)
  out <- data.frame(Term = unname(lab), row.names = NULL)
  for (m in models) {
    x <- coefs[coefs$model == m, ]
    cell <- setNames(sprintf("%s [%s, %s]%s", fmt(x$OR), fmt(x$OR_lo), fmt(x$OR_hi), stars(x$p)), x$term)
    out[[m]] <- ifelse(terms %in% names(cell), cell[terms], "")
  }
  out <- out[apply(out[, -1, drop = FALSE], 1, function(r) any(r != "")), ]
  for (r in c("people", "person_years", "events"))
    out <- rbind(out, c(r, sapply(models, function(m) coefs[[r]][coefs$model == m][1])))
  out
}
tA <- mk_table(c("Window 1: Wave I + adversity",
                 "Window 1: no self-reported juvenile delinquency",
                 "Window 1: no SERIOUS juvenile delinquency"))
tB <- mk_table(c("Window 2: Wave I + adversity", "Window 2: + Wave III roles",
                 "Window 2: + role quality", "Window 2: no self-reported juvenile delinquency",
                 "Window 2: no SERIOUS juvenile delinquency"))
tF <- mk_table("Pooled windows: predictor x window")
tG <- mk_table(c("Reported ages only: 18 to Wave IV", "Reported ages: + ADHD x age",
                 "Reported ages: + ADHD x age + peer x age"))
split_names <- function(tag) paste(tag, c("split: perceived peer substance use", "split: + own group offending",
                                            "split: friends' own reports", "split: all four"))
tI <- mk_table(split_names("Window 1")); tJ <- mk_table(split_names("Window 2"))
tH <- mk_table(c("Window 1: Wave I + adversity", "Window 1: + own Wave I delinquency variety",
                 "Window 1: self-reported non-offenders only",
                 "Window 2: + Wave III roles", "Window 2: + own Wave III offending"))

# Where the failures fall, and how many had undetected juvenile offending
fails <- d[d$in_win1 & !is.na(d$age_first), ]
tC <- data.frame(
  Group = c("First arrest 18-24", "First arrest 25+", "  of which age imputed (W4-W5)",
            "BROAD screen: onset with NO Wave I/II self-reported delinquency",
            "BROAD screen: onset WITH Wave I/II self-reported delinquency",
            "STRICT screen: onset with NO serious Wave I/II delinquency",
            "STRICT screen: onset WITH serious Wave I/II delinquency",
            "Juvenile self-report unknown (broad screen)"),
  n = c(sum(fails$age_first < 25), sum(fails$age_first >= 25), sum(fails$age_first_imputed),
        sum(fails$sr_juv_delinq %in% 0), sum(fails$sr_juv_delinq %in% 1),
        sum(fails$sr_juv_strict %in% 0), sum(fails$sr_juv_strict %in% 1),
        sum(is.na(fails$sr_juv_delinq))))
tC$pct_of_onset <- fmt(100 * tC$n / sum(fails$age_first >= 18), 1)
# Same split among people who never had an adult arrest, for comparison
noarr <- d[d$in_win1 & is.na(d$age_first), ]
tC$comparison_no_arrest_pct <- c("--", "--", "--",
  fmt(100 * mean(noarr$sr_juv_delinq %in% 0), 1), fmt(100 * mean(noarr$sr_juv_delinq %in% 1), 1),
  fmt(100 * mean(noarr$sr_juv_strict %in% 0), 1), fmt(100 * mean(noarr$sr_juv_strict %in% 1), 1), "--")

tD <- data.frame(Measure = lab[w3_block],
  pct_at_w3 = sapply(w3_block, function(v) fmt(100 * weighted.mean(d[[v]][d$in_win2], d$GSW1345[d$in_win2]), 1)),
  n_yes = sapply(w3_block, function(v) sum(d[[v]][d$in_win2] == 1)), row.names = NULL)

tabs <- list("Table A. Window 1 (expiry): first arrest between 18 and the Wave III interview" = tA,
             "Table B. Window 2 (late): first arrest after the Wave III interview" = tB,
             "Table C. Where the adult-onset arrests fall" = tC,
             "Table D. Wave III roles in the window 2 risk set (weighted %)" = tD,
             "Table E. Pooled windows: does early adversity act differently after Wave III?" = tF,
             "Table F. Reported ages only (no imputed arrest ages): age gradient and ADHD x age" = tG,
             "Table G. Detection vs. involvement: peer deviance net of own reported offending" = tH,
             "Table H. Splitting peer deviance, window 1 (network subsample)" = tI,
             "Table I. Splitting peer deviance, window 2 (network subsample)" = tJ,
             "Table J. Sample flow" = flow_log)
fn <- c("tableA_window1", "tableB_window2", "tableC_onset_composition", "tableD_w3_roles",
        "tableE_pooled_window_test", "tableF_reported_ages", "tableG_detection_vs_involvement",
        "tableH_split_window1", "tableI_split_window2", "tableJ_flow")
for (i in seq_along(tabs)) write.csv(tabs[[i]], file.path(ao_dir, paste0(fn[i], ".csv")), row.names = FALSE)
ht <- function(df) paste0("<table><tr>", paste0("<th>", names(df), "</th>", collapse = ""), "</tr>",
  paste(apply(df, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>")), collapse = ""), "</table>")
writeLines(c("<html><head><meta charset='utf-8'><style>body{font-family:Georgia,serif;margin:30px;max-width:1300px}",
  "table{border-collapse:collapse;font-size:13px;margin-bottom:28px}th,td{border:1px solid #bbb;padding:4px 8px}",
  "th{background:#eee}</style></head><body><h1>Adult-onset arrest: two stacked risk sets</h1>",
  "<p>Discrete-time logistic hazard models on person-years; odds ratios [95% CI]; * p&lt;.05 ** p&lt;.01 *** p&lt;.001.",
  "Weighted (GSW1345) with CLUSTER2 cluster-robust standard errors.</p>",
  unlist(lapply(names(tabs), function(n) c(paste0("<h2>", n, "</h2>"), ht(tabs[[n]])))),
  "</body></html>"), file.path(ao_dir, "ARREST_ONSET.html"))
print(tC, row.names = FALSE)

# ---------------------------------------------------------------------------
# 8. FIGURES
# ---------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 11, base_family = "serif") +
          theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
                legend.position = "bottom"))
# Age-specific hazard of first arrest, REPORTED AGES ONLY. Ages with thin
# exposure are cut, because a handful of person-years produces noise that
# looks like structure.
hz <- aggregate(cbind(e = event * GSW1345, w = GSW1345) ~ age, data = pyr, sum)
hz$py <- as.vector(table(pyr$age)[as.character(hz$age)])
hz$hazard <- 100 * hz$e / hz$w
hz <- hz[hz$py >= 500, ]
p1 <- ggplot(hz, aes(age, hazard)) + geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Age", y = "First arrests per 100 still at risk",
       title = "Figure A. Age-specific hazard of first adult arrest",
       subtitle = "Reported ages only; everyone censored at the Wave IV interview",
       caption = paste("Weighted (GSW1345). The 229 Wave IV-to-Wave V arrests with no reported age are",
         "excluded rather than\nimputed. Ages with fewer than 500 person-years are not plotted.",
         "The decline is also consistent with\nselection (higher-risk people fail out of the risk set first),",
         "not only with maturation."))
ggsave(file.path(ao_dir, "figA_hazard_by_age.png"), p1, width = 7.5, height = 4.6, dpi = 300, bg = "white")

plot_or <- function(models, terms, title, file, w = 8) {
  f <- coefs[coefs$model %in% models & coefs$term %in% terms, ]
  f$term_lab <- factor(lab[f$term], levels = rev(lab[terms]))
  f$model <- factor(f$model, levels = models)
  p <- ggplot(f, aes(OR, term_lab)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(xmin = OR_lo, xmax = OR_hi), size = 0.3) +
    facet_wrap(~ model, nrow = 1) + scale_x_log10() +
    labs(x = "Odds ratio (log scale), 95% CI", y = NULL, title = title)
  ggsave(file.path(ao_dir, file), p, width = w, height = 4.4, dpi = 300, bg = "white")
}
plot_or(c("Window 1: Wave I + adversity", "Window 2: + Wave III roles"),
        c("strain", "peer", "bond", "parent_control_z", "phys_abuse", "sex_abuse", "adhd"),
        "Figure B. Adolescent strain, bonds and adversity, by window", "figB_windows.png")
plot_or(c("Window 2: + Wave III roles", "Window 2: + role quality"),
        c(w3_block, "working_w3:hardship_w3", "partner_w3:hardship_w3"),
        "Figure C. Adult roles at Wave III and their quality", "figC_roles.png")

saveRDS(list(d = d, coefs = coefs, lab = lab, py1 = py1, py2 = py2, w3_block = w3_block),
        file.path(work_dir, "arrest_onset_results.rds"))
cat("\nSaved tables and figures to output/arrest_onset/\n")

# ---------------------------------------------------------------------------
# 9. HEADLINE RESULTS -- printed so the run log records the findings, not
#    just the sample sizes. Everything here is also in the tables.
# ---------------------------------------------------------------------------
show <- function(model, term, note = "") {
  # R writes interactions in whichever order the formula produced, so try both
  alt <- if (grepl(":", term)) paste(rev(strsplit(term, ":")[[1]]), collapse = ":") else term
  x <- coefs[coefs$model == model & coefs$term %in% c(term, alt), ]
  if (!nrow(x)) return(invisible(NULL))
  cat(sprintf("  %-46s %-26s OR %5.2f [%4.2f, %5.2f] %-3s  (events %d)%s\n",
              substr(model, 1, 46), lab[[term]], x$OR, x$OR_lo, x$OR_hi,
              stars(x$p), x$events, note))
}
cat("\n================ STAGE 5: HEADLINE RESULTS ================\n")
cat("\n-- Age gradient and the maturation test (reported ages only) --\n")
show("Reported ages: + ADHD x age", "adhd:age_c", "  <- flat = age decline does not track self-regulation")
show("Reported ages only: 18 to Wave IV", "age_c")
show("Reported ages only: 18 to Wave IV", "adhd")

cat("\n-- Does early adversity act differently after Wave III? --\n")
for (t in c("late:sex_abuse", "late:phys_abuse", "late:adhd", "late:peer", "late:strain"))
  show("Pooled windows: predictor x window", t)

cat("\n-- Peer deviance across every screen (the one that holds) --\n")
for (m in c("Window 1: Wave I + adversity", "Window 1: no self-reported juvenile delinquency",
            "Window 1: no SERIOUS juvenile delinquency", "Window 2: + Wave III roles",
            "Window 2: no self-reported juvenile delinquency",
            "Window 2: no SERIOUS juvenile delinquency",
            "Reported ages only: 18 to Wave IV")) show(m, "peer")

cat("\n-- Detection vs involvement: peer deviance net of own offending --\n")
show("Window 1: Wave I + adversity", "peer", "  <- baseline")
show("Window 1: + own Wave I delinquency variety", "peer")
show("Window 1: + own Wave I delinquency variety", "delinq_variety_w1")
show("Window 1: self-reported non-offenders only", "peer", "  <- no self-reported juvenile offending")
show("Window 2: + Wave III roles", "peer", "  <- baseline")
show("Window 2: + own Wave III offending", "peer")
show("Window 2: + own Wave III offending", "offend_w3_c")

cat("\n-- Splitting peer deviance: substance use vs group offending (network subsample) --\n")
for (tag in c("Window 1", "Window 2")) {
  m <- split_names(tag)
  show(m[1], "peer", "  <- baseline on this subsample")
  show(m[2], "peer"); show(m[2], "group_fight")
  show(m[3], "fr_su_z"); show(m[3], "fr_risk_z")
  show(m[4], "peer"); show(m[4], "fr_su_z"); show(m[4], "fr_risk_z"); show(m[4], "group_fight")
  cat("\n")
}

cat("\n-- Adult roles at Wave III (window 2) and their quality --\n")
for (t in c(w3_block, "working_w3:hardship_w3", "partner_w3:hardship_w3"))
  show(if (grepl(":", t)) "Window 2: + role quality" else "Window 2: + Wave III roles", t)

cat("\n-- How much of official adult onset is undetected juvenile offending --\n")
cat(sprintf("  BROAD screen:  %s%% of onset cases reported juvenile delinquency (vs %s%% of the never-arrested)\n",
            tC$pct_of_onset[5], tC$comparison_no_arrest_pct[5]))
cat(sprintf("  STRICT screen: %s%% reported SERIOUS juvenile delinquency (vs %s%% of the never-arrested)\n",
            tC$pct_of_onset[7], tC$comparison_no_arrest_pct[7]))
cat(sprintf("  %d of %d adult-onset arrests (%s%%) have an imputed age (Wave IV-V, no age recorded)\n",
            tC$n[3], sum(tC$n[1:2]), tC$pct_of_onset[3]))
cat("\nFull detail: output/arrest_onset/ARREST_ONSET.html\n")
