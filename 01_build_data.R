###############################################################################
# 01_build_data.R
# STAGE 1 -- BUILD THE ANALYSIS DATA
#
# Reads only the needed columns from the raw ICPSR files (Waves I, III, IV,
# V and weights), builds every variable defined in 00_setup.R, and prints the
# sample size after EVERY filtering step.
#
# Produces:
#   data_work/analysis_data.rds          one row per Wave I respondent
#   output/logs/sample_flow.csv          the N-by-step table
###############################################################################

paths <- locate_icpsr()
cat("\n================ STAGE 1: DATA BUILD ================\n")

# Helper for "marked / not marked" charge items (7 = legitimate skip -> 0)
chg <- function(x, skip_is_zero = TRUE) {
  x <- clean_num(x, c(6, 8, 9))
  if (skip_is_zero) x[x == 7] <- 0
  x
}
any_of_cols <- function(m) as.numeric(rowSums(m == 1, na.rm = TRUE) > 0)

# ---------------------------------------------------------------------------
# 1. WAVE I
# ---------------------------------------------------------------------------
su_items <- c("H1TO17", "H1TO30", "H1TO34", "H1TO37", "H1TO40")
w1_vars <- unique(c("AID", "IMONTH", "IYEAR", "H1GI1M", "H1GI1Y", "BIO_SEX",
                    "H1GI4", paste0("H1GI6", LETTERS[1:5]), "H1RM1", "H1RF1",
                    names(w1_delinq_items), cesd_items, peer_items,
                    school_items, family_items, su_items))
w1 <- as.data.frame(read_dta(paths$w1_main, col_select = all_of(w1_vars)))
cat("Wave I in-home file read:", nrow(w1), "respondents\n")

# --- Wave I delinquency (0/1) ------------------------------------------------
for (v in names(w1_delinq_items)) w1[[v]] <- clean_num(w1[[v]], c(6, 8, 9))
w1$delinq_w1 <- any_act(w1, names(w1_delinq_items))
w1$delinq_variety_w1 <- rowSums(w1[, names(w1_delinq_items)] > 0, na.rm = TRUE)
w1$delinq_variety_w1[is.na(w1$delinq_w1)] <- NA

# --- Strain: CES-D ------------------------------------------------------------
for (v in cesd_items) {
  w1[[v]] <- clean_num(w1[[v]], c(6, 8, 9))
  if (v %in% cesd_reverse) w1[[v]] <- 3 - w1[[v]]
}
n_cesd <- rowSums(!is.na(w1[, cesd_items]))
w1$strain_raw <- ifelse(n_cesd >= 15,
                        rowMeans(w1[, cesd_items], na.rm = TRUE) * 19, NA)

# --- Peer deviance ------------------------------------------------------------
for (v in peer_items) w1[[v]] <- clean_num(w1[[v]], c(6, 7, 8, 9))
w1$peer_raw <- rowSums(w1[, peer_items])

# --- Bond deficit (higher = weaker bonds) -----------------------------------
for (v in school_items) w1[[v]] <- clean_num(w1[[v]], c(6, 7, 8))
for (v in family_items) {
  w1[[v]] <- clean_num(w1[[v]], c(6, 96, 98))
  w1[[v]] <- 6 - w1[[v]]
}
bond_items <- c(school_items, family_items)
bond_z <- sapply(bond_items, function(v) zscore(w1[[v]]))
w1$bond_raw <- ifelse(rowSums(!is.na(bond_z)) >= 4,
                      rowMeans(bond_z, na.rm = TRUE), NA)

# --- Own Wave I substance use (robustness Model D) --------------------------
binge <- clean_num(w1$H1TO17, c(96, 98, 99))       # 1-6 = any binge; 7 never; 97 never drank
w1$w1_binge <- ifelse(binge %in% 1:6, 1, ifelse(binge %in% c(7, 97), 0, NA))
agefirst <- function(v) { x <- clean_num(v, c(96, 98, 99)); ifelse(x > 0, 1, x) }
w1$w1_marij <- agefirst(w1$H1TO30)
w1$w1_hard  <- pmax(agefirst(w1$H1TO34), agefirst(w1$H1TO37), agefirst(w1$H1TO40))
w1$w1_su_count <- w1$w1_binge + w1$w1_marij + w1$w1_hard

# --- Controls -----------------------------------------------------------------
bm <- clean_num(w1$H1GI1M, 96); by <- clean_num(w1$H1GI1Y, 96)
w1$age_w1 <- ((clean_num(w1$IYEAR, NA) * 12 + clean_num(w1$IMONTH, NA)) -
              (by * 12 + bm)) / 12
w1$male <- ifelse(clean_num(w1$BIO_SEX, 6) == 1, 1, 0)
hisp_raw <- clean_num(w1$H1GI4, c(6, 8))
race <- sapply(paste0("H1GI6", LETTERS[1:5]), function(v) clean_num(w1[[v]], c(6, 8)))
w1$race_eth <- ifelse(hisp_raw == 1, "Hispanic",
               ifelse(race[, "H1GI6B"] == 1, "NH Black",
               ifelse(race[, "H1GI6A"] == 1 & rowSums(race[, -1], na.rm = TRUE) == 0,
                      "NH White", "NH Other")))
w1$race_eth[is.na(hisp_raw) | rowSums(is.na(race)) == 5] <- NA
w1$hisp       <- as.numeric(w1$race_eth == "Hispanic")
w1$black      <- as.numeric(w1$race_eth == "NH Black")
w1$other_race <- as.numeric(w1$race_eth == "NH Other")
pe <- cbind(clean_num(w1$H1RM1, c(11, 12, 96, 97, 98, 99)),
            clean_num(w1$H1RF1, c(11, 12, 96, 97, 98, 99)))
w1$par_college <- ifelse(rowSums(pe >= 8 & pe <= 9, na.rm = TRUE) > 0, 1, 0)
w1$par_ed_unknown <- ifelse(w1$par_college == 0 & rowSums(!is.na(pe)) == 0, 1, 0)

# Standardized focal predictors (standardized on all valid Wave I cases)
w1$strain <- zscore(w1$strain_raw)
w1$peer   <- zscore(w1$peer_raw)
w1$bond   <- zscore(w1$bond_raw)
w1$w1_su  <- zscore(w1$w1_su_count)

w1 <- w1[, c("AID", "delinq_w1", "delinq_variety_w1", "strain_raw", "peer_raw",
             "bond_raw", "strain", "peer", "bond", "w1_binge", "w1_marij",
             "w1_hard", "w1_su_count", "w1_su", "age_w1", "male", "race_eth",
             controls[!controls %in% c("age_w1", "male")])]

# ---------------------------------------------------------------------------
# 2. WEIGHTS
# ---------------------------------------------------------------------------
w1w <- as.data.frame(read_dta(paths$w1_weight))
w5w <- as.data.frame(read_sas(paths$w5_weight, col_select = c("AID", "GSW145", "GSW1345")))

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 3. WAVE III: offending, illicit drug use, childhood ADHD symptoms
# ---------------------------------------------------------------------------
w3_drug_since95 <- c("H3TO108", "H3TO111", "H3TO114", "H3TO117", "H3TO120")   # any use since 6/1995
w3_rx_since95   <- paste0("H3TO105", LETTERS[1:4])                           # nonmedical Rx since 6/1995
w3_drug_pastyr  <- c(H3TO109 = "H3TO108", H3TO112 = "H3TO111",               # past-year item = gate item
                     H3TO115 = "H3TO114", H3TO118 = "H3TO117")
w3 <- as.data.frame(read_dta(paths$w3_main, col_select = all_of(c(
  "AID", "CALCAGE3", names(w3_offend_items), w3_drug_since95, w3_rx_since95,
  names(w3_drug_pastyr), adhd_inatt, adhd_hyper))))
cat("Wave III in-home file read:", nrow(w3), "respondents\n")

for (v in setdiff(names(w3_offend_items), "H3DS17")) w3[[v]] <- clean_num(w3[[v]], c(6, 8, 9))
w3$H3DS17 <- clean_num(w3$H3DS17, c(996, 998, 999))
w3$offend_w3 <- any_act(w3, names(w3_offend_items))
w3$age_w3 <- clean_num(w3$CALCAGE3, numeric(0))

# Drug use since June 1995 (used for the late-onset screen)
since <- sapply(c(w3_drug_since95, w3_rx_since95), function(v) {
  x <- clean_num(w3[[v]], c(6, 8, 9)); x[x == 7] <- 0; x })
w3$drug_since95_w3 <- any_act(as.data.frame(since), colnames(since))
# Past-year illicit drug use (used as an onset window in the secondary design)
pastyr <- sapply(names(w3_drug_pastyr), function(v) {
  x <- clean_num(w3[[v]], c(6, 8, 9))
  gate <- clean_num(w3[[w3_drug_pastyr[[v]]]], c(6, 8, 9))
  x[x == 7 & gate %in% 0] <- 0; x[x == 7] <- NA; x })
w3$drug_pastyr_w3 <- any_act(as.data.frame(pastyr), colnames(pastyr))
# Same, excluding marijuana (H3TO109), for the no-marijuana sensitivity
w3$hard_pastyr_w3 <- any_act(as.data.frame(pastyr), setdiff(colnames(pastyr), "H3TO109"))

# Childhood ADHD symptoms (ages 5-12, recalled at Wave III)
adhd_m <- sapply(c(adhd_inatt, adhd_hyper), function(v) clean_num(w3[[v]], c(6, 8)))
often <- adhd_m >= 2
n_ans <- rowSums(!is.na(adhd_m))
w3$adhd_count <- ifelse(n_ans >= 14, rowSums(often, na.rm = TRUE) * 17 / n_ans, NA)
w3$adhd_probable <- ifelse(n_ans >= 14,
  as.numeric(rowSums(often[, adhd_inatt], na.rm = TRUE) >= 6 |
             rowSums(often[, adhd_hyper], na.rm = TRUE) >= 6), NA)
w3$in_w3 <- 1
w3 <- w3[, c("AID", "in_w3", "age_w3", "offend_w3", "drug_since95_w3",
             "drug_pastyr_w3", "hard_pastyr_w3", "adhd_count", "adhd_probable")]

# ---------------------------------------------------------------------------
# 4. WAVE IV: offending, illicit drug use, age at first use, childhood abuse
# ---------------------------------------------------------------------------
w4 <- as.data.frame(read_dta(paths$w4_main, col_select = all_of(c(
  "AID", names(w4_offend_items), "H4TO65B", "H4TO68", "H4TO70", "H4TO96", "H4TO98",
  "H4MA3", "H4MA5"))))
cat("Wave IV in-home file read:", nrow(w4), "respondents\n")

for (v in names(w4_offend_items)) w4[[v]] <- clean_num(w4[[v]], c(6, 8))
w4$H4DS12[w4$H4DS12 == 7] <- 0
w4$offend_w4 <- any_act(w4, names(w4_offend_items))

# Past-year use: marijuana days (H4TO70) and most-used other illegal drug
# days (H4TO98). 97 = legitimate skip = never used that drug -> 0.
mj_days  <- clean_num(w4$H4TO70, c(96, 98)); mj_days[mj_days == 97] <- 0
oth_days <- clean_num(w4$H4TO98, c(96, 98)); oth_days[oth_days == 97] <- 0
w4$mj_pastyr_w4   <- as.numeric(mj_days >= 1)
w4$hard_pastyr_w4 <- as.numeric(oth_days >= 1)
w4$drug_pastyr_w4 <- ifelse(w4$mj_pastyr_w4 %in% 1 | w4$hard_pastyr_w4 %in% 1, 1,
                     ifelse(!is.na(w4$mj_pastyr_w4) & !is.na(w4$hard_pastyr_w4), 0, NA))
w4$ever_mj_w4 <- clean_num(w4$H4TO65B, c(6, 8))
# Juvenile drug use screen: first marijuana use, or first use of the
# most-used other drug, before age 18 (97 = never used -> not juvenile).
age_mj  <- clean_num(w4$H4TO68, c(96, 98))
age_oth <- clean_num(w4$H4TO96, c(96, 98))
juv <- cbind(ifelse(age_mj == 97, 0, as.numeric(age_mj < 18)),
             ifelse(age_oth == 97, 0, as.numeric(age_oth < 18)))
w4$juv_drug_w4 <- ifelse(rowSums(juv == 1, na.rm = TRUE) > 0, 1,
                  ifelse(rowSums(is.na(juv)) == 0, 0, NA))

abuse <- function(x) { x <- clean_num(x, c(96, 98)); ifelse(x == 6, 0, ifelse(x %in% 1:5, 1, NA)) }
w4$phys_abuse <- abuse(w4$H4MA3)
w4$sex_abuse  <- abuse(w4$H4MA5)
w4$in_w4 <- 1
w4 <- w4[, c("AID", "in_w4", "offend_w4", "mj_pastyr_w4", "hard_pastyr_w4",
             "drug_pastyr_w4", "ever_mj_w4", "juv_drug_w4", "phys_abuse", "sex_abuse")]

# ---------------------------------------------------------------------------
# 5. WAVE V: offending, illicit drug use
# ---------------------------------------------------------------------------
w5_hard <- c("H5TO27A", "H5TO27B", "H5TO27C", "H5TO27D")   # past-30-day cocaine, meth, heroin, other
w5_rx   <- paste0("H5TO26", LETTERS[1:4])                   # nonmedical prescription drugs
w5 <- as.data.frame(read_dta(paths$w5_main, col_select = all_of(c(
  "AID", names(w5_offend_items), "H5TO20", "H5TO21", w5_hard, w5_rx))))
cat("Wave V survey file read:", nrow(w5), "respondents\n")
for (v in names(w5_offend_items)) w5[[v]] <- clean_num(w5[[v]], numeric(0))
w5$offend_w5 <- any_act(w5, names(w5_offend_items))
mj30 <- clean_num(w5$H5TO21, numeric(0)); mj30[mj30 == 997] <- 0; mj30[mj30 == 97] <- 0
w5$ever_mj_w5 <- clean_num(w5$H5TO20, numeric(0))
w5$mj_30d_w5 <- as.numeric(mj30 >= 1)
w5$hard_w5 <- any_act(w5, c(w5_hard, w5_rx))
w5$in_w5 <- 1
w5 <- w5[, c("AID", "in_w5", "offend_w5", "ever_mj_w5", "mj_30d_w5", "hard_w5")]

# ---------------------------------------------------------------------------
# 6. MERGE (Wave I is the base: every Wave I respondent keeps a row)
# ---------------------------------------------------------------------------
dat <- merge(w1, w1w, by = "AID", all.x = TRUE)
dat <- merge(dat, w5w, by = "AID", all.x = TRUE)
dat <- merge(dat, w3,  by = "AID", all.x = TRUE)
dat <- merge(dat, w4,  by = "AID", all.x = TRUE)
dat <- merge(dat, w5,  by = "AID", all.x = TRUE)
for (v in c("in_w3", "in_w4", "in_w5")) dat[[v]][is.na(dat[[v]])] <- 0
dat$adhd <- zscore(dat$adhd_count)

# Wave I illicit drug use (ever marijuana, cocaine, inhalants, other drugs)
dat$w1_drug <- ifelse(dat$w1_marij %in% 1 | dat$w1_hard %in% 1, 1,
               ifelse(!is.na(dat$w1_marij) & !is.na(dat$w1_hard), 0, NA))
# First-ever marijuana use between Waves IV and V (adult by definition)
dat$new_mj_w5 <- ifelse(dat$ever_mj_w4 == 1, 0,
                 ifelse(dat$ever_mj_w4 == 0, dat$ever_mj_w5, NA))

# ---------------------------------------------------------------------------
# 7. ONSET OUTCOMES
# ---------------------------------------------------------------------------
onset_combine <- function(parts) {
  m <- as.matrix(parts)
  ifelse(rowSums(m == 1, na.rm = TRUE) > 0, 1, ifelse(rowSums(is.na(m)) == 0, 0, NA))
}
adult_w3 <- function(x) ifelse(dat$age_w3 >= 19, x, NA)   # W3 counts only if 19+

# Offending (self-report only)
dat$off_late <- onset_combine(dat[, c("offend_w4", "offend_w5")])
dat$off_all  <- onset_combine(cbind(adult_w3(dat$offend_w3), dat$offend_w4, dat$offend_w5))
# Illicit drug use (any, incl. marijuana)
dat$sub_late <- onset_combine(dat[, c("drug_pastyr_w4", "mj_30d_w5", "new_mj_w5", "hard_w5")])
dat$sub_all  <- onset_combine(cbind(adult_w3(dat$drug_pastyr_w3), dat$drug_pastyr_w4,
                                    dat$mj_30d_w5, dat$new_mj_w5, dat$hard_w5))
# Illicit drug use EXCLUDING marijuana (sensitivity: marijuana was legalized
# for adults in several states during the Wave V field period)
dat$hard_late <- onset_combine(dat[, c("hard_pastyr_w4", "hard_w5")])
dat$hard_all  <- onset_combine(cbind(adult_w3(dat$hard_pastyr_w3), dat$hard_pastyr_w4, dat$hard_w5))
# Any adult onset = offending OR illicit drug use
dat$onset_late <- onset_combine(dat[, c("off_late", "sub_late")])
dat$onset_all  <- onset_combine(dat[, c("off_all", "sub_all")])
dat$onset_w45  <- dat$onset_late                     # sensitivity design uses the W4-W5 window

route <- function(off, sub) ifelse(is.na(off) | is.na(sub), NA,
  ifelse(off == 1 & sub == 1, "both", ifelse(off == 1, "offending_only",
  ifelse(sub == 1, "substance_only", "none"))))
dat$route_late <- route(dat$off_late, dat$sub_late)
dat$route_all  <- route(dat$off_all,  dat$sub_all)

# ---------------------------------------------------------------------------
# 8. SAMPLE FLOW -- printed and saved at every step
# ---------------------------------------------------------------------------
complete_x   <- complete.cases(dat[, c("strain", "peer", "bond", controls)])
complete_adv <- complete.cases(dat[, adversity])

cat("\n---- Sample flow: H1 (Wave I delinquency) ----\n")
s <- rep(TRUE, nrow(dat))
record_step("H1", "Wave I public-use in-home respondents", sum(s))
s <- s & !is.na(dat$GSWGT1);        record_step("H1", "has Wave I weight (GSWGT1)", sum(s))
s <- s & !is.na(dat$delinq_w1);     record_step("H1", "valid Wave I delinquency", sum(s))
s <- s & !is.na(dat$strain);        record_step("H1", "valid strain (CES-D)", sum(s))
s <- s & !is.na(dat$peer);          record_step("H1", "valid peer deviance", sum(s))
s <- s & !is.na(dat$bond);          record_step("H1", "valid bond deficit", sum(s))
s <- s & complete.cases(dat[, controls])
record_step("H1", "valid controls = H1 ANALYTIC SAMPLE", sum(s))
dat$in_h1 <- s

onset_flow <- function(label, outcome, wt, use_w3) {
  cat("\n---- Sample flow:", label, "----\n")
  s <- rep(TRUE, nrow(dat))
  record_step(label, "Wave I public-use in-home respondents", sum(s))
  s <- s & dat$delinq_w1 %in% 0;    record_step(label, "Wave I: no delinquent acts", sum(s))
  s <- s & dat$w1_drug %in% 0;      record_step(label, "Wave I: never used illicit drugs", sum(s))
  if (use_w3 != "none") { s <- s & dat$in_w3 == 1; record_step(label, "interviewed at Wave III", sum(s)) }
  s <- s & dat$in_w4 == 1;          record_step(label, "interviewed at Wave IV", sum(s))
  s <- s & dat$in_w5 == 1;          record_step(label, "interviewed at Wave V", sum(s))
  s <- s & !is.na(dat[[wt]]);       record_step(label, paste0("has longitudinal weight (", wt, ")"), sum(s))
  s <- s & dat$juv_drug_w4 %in% 0
  record_step(label, "first drug use not before 18 (Wave IV report)", sum(s))
  if (use_w3 == "screen") {
    s <- s & dat$offend_w3 %in% 0 & dat$drug_since95_w3 %in% 0
    record_step(label, "no offending or drug use by Wave III (late onset)", sum(s))
  }
  if (use_w3 == "window") {
    s <- s & !((dat$offend_w3 %in% 1 | dat$drug_pastyr_w3 %in% 1) & dat$age_w3 < 19)
    record_step(label, "drop Wave III offenders/users aged 18 (partly juvenile)", sum(s))
  }
  s <- s & !is.na(dat[[outcome]]);  record_step(label, "valid onset outcome", sum(s))
  s <- s & complete_x;              record_step(label, "valid predictors + controls", sum(s))
  if (use_w3 != "none") {
    s <- s & complete_adv
    record_step(label, "valid abuse + ADHD", sum(s))
  }
  s <- s & !is.na(dat$w1_binge);    record_step(label, "valid W1 binge drinking = ANALYTIC SAMPLE", sum(s))
  s
}
# PRIMARY uses Wave III as an onset window (all adult onset, ~ages 19-43);
# SECONDARY keeps Wave III as a screen (late onset only, ~ages 25-43).
dat$in_h2p <- onset_flow("H2 primary (all adult onset)",    "onset_all",  "GSW1345", "window")
dat$in_h2s <- onset_flow("H2 secondary (late onset only)",  "onset_late", "GSW1345", "screen")
dat$in_h2x <- onset_flow("H2 sensitivity (W1 screen only)", "onset_w45",  "GSW145",  "none")

cat("\nOnset in analytic samples (unweighted):\n")
for (f in list(c("in_h2p", "onset_all", "route_all"), c("in_h2s", "onset_late", "route_late"),
               c("in_h2x", "onset_w45", "route_late"))) {
  y <- dat[[f[2]]][dat[[f[1]]]]
  cat(sprintf("  %-7s %d of %d (%.1f%%)\n", f[1], sum(y), length(y), 100 * mean(y)))
  print(table(dat[[f[3]]][dat[[f[1]]]]))
}

# 9. SAVE
# ---------------------------------------------------------------------------
saveRDS(dat, file.path(work_dir, "analysis_data.rds"))
write.csv(flow_log, file.path(log_dir, "sample_flow.csv"), row.names = FALSE)
write.csv(weight_plan, file.path(log_dir, "weight_plan.csv"), row.names = FALSE)
cat("\nSaved data_work/analysis_data.rds and output/logs/sample_flow.csv\n")

rm(w1, w1w, w3, w4, w5, w5w, bond_z, race, pe, adhd_m, often, since, pastyr, juv); invisible(gc())
