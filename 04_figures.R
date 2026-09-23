###############################################################################
# 04_figures.R
# STAGE 4 -- FIGURES (all drawn straight from the data and fitted models)
#
# Figure 1  Sample flow: N after each filtering step, every sample
# Figure 2  Wave I delinquency by strain at low/average/high peer deviance
#           (H1, weighted Model B)
# Figure 3  Adult onset by strain at low/average/high peer deviance
#           (H2 primary, weighted Model C)
# Figure 4  Adult onset: offending only, drug use only, both (weighted %)
# Figure 5  Weighted vs. unweighted odds ratios, focal and adversity terms
# Figure 6  Offending onset vs. illicit drug use onset: odds ratios (Model C)
#
# Each figure: 300-dpi PNG + vector PDF in output/figures/.
###############################################################################

cat("\n================ STAGE 4: FIGURES ================\n")
dat  <- readRDS(file.path(work_dir, "analysis_data.rds"))
mr   <- readRDS(file.path(work_dir, "model_results.rds"))
flow <- read.csv(file.path(log_dir, "sample_flow.csv"))
fits <- mr$fits; ct <- mr$coef_table

theme_set(theme_minimal(base_size = 11, base_family = "serif") +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold"),
                legend.position = "bottom"))
save_fig <- function(p, name, w = 7, h = 4.5) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  cat("  saved", name, "\n")
}
unlink(list.files(fig_dir, full.names = TRUE))          # clear figures from earlier runs

# ---------------------------------------------------------------------------
# FIGURE 1: SAMPLE FLOW
# ---------------------------------------------------------------------------
flow$step_no <- ave(seq_len(nrow(flow)), flow$sample, FUN = seq_along)
flow$lab <- paste0(flow$step_no, ". ", flow$step)
flow$sample <- factor(flow$sample, levels = unique(flow$sample))
flow$key <- paste(flow$sample, sprintf("%02d", flow$step_no))
lab_map <- setNames(flow$lab, flow$key)
p1 <- ggplot(flow, aes(x = n, y = reorder(key, -step_no))) +
  geom_col(fill = "grey55", width = 0.7) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.1, size = 2.8, family = "serif") +
  facet_wrap(~ sample, ncol = 1, scales = "free_y") +
  scale_y_discrete(labels = function(k) lab_map[k]) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Respondents remaining", y = NULL,
       title = "Figure 1. Sample flow from the Wave I public-use file")
save_fig(p1, "figure1_sample_flow", w = 8.5, h = 13)

# ---------------------------------------------------------------------------
# FIGURES 2 AND 3: PREDICTED PROBABILITIES (strain x peer deviance)
# Other covariates held at weighted sample means; ribbons = 95% CI.
# ---------------------------------------------------------------------------
pred_curves <- function(fit, d, wt, hold) {
  grid <- expand.grid(strain = seq(-2, 2, by = 0.1), peer = c(-1, 0, 1))
  for (v in hold) grid[[v]] <- sum(d[[v]] * d[[wt]]) / sum(d[[wt]])
  pr <- predict(fit, newdata = grid, type = "link", se.fit = TRUE)
  eta <- as.numeric(coef(pr)); se <- as.numeric(SE(pr))
  grid$prob <- plogis(eta); grid$lo <- plogis(eta - 1.96 * se); grid$hi <- plogis(eta + 1.96 * se)
  grid$peer_lab <- factor(grid$peer, levels = c(-1, 0, 1),
    labels = c("Low peer deviance (-1 SD)", "Average (mean)", "High peer deviance (+1 SD)"))
  grid
}
plot_curves <- function(g, ylab, title, caption) {
  ggplot(g, aes(strain, prob, linetype = peer_lab, fill = peer_lab)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) + scale_fill_grey(start = 0.2, end = 0.6) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    labs(x = "Strain (CES-D, standard deviations from mean)", y = ylab,
         linetype = NULL, fill = NULL, title = title, caption = caption)
}
int_caption <- function(spec, model, extra) {
  b <- ct[ct$spec == spec & ct$model == model & ct$weighting == "weighted" &
          ct$family == "logit" & ct$term == "strain:peer", ]
  sprintf("%s\nInteraction b = %.2f (SE %.2f), p = %.3f. N = %d.", extra, b$b, b$se, b$p, b$n)
}

h1 <- dat[dat$in_h1, ]
g2 <- pred_curves(fits$H1_B_weighted, h1, "GSWGT1", c("bond", controls))
save_fig(plot_curves(g2, "Predicted probability of delinquency",
  "Figure 2. Strain and peer deviance: Wave I delinquency",
  int_caption("H1", "B", "Weighted logistic Model B (GSWGT1, CLUSTER2).")), "figure2_H1_interaction")

h2p <- dat[dat$in_h2p, ]
g3 <- pred_curves(fits$H2P_C_weighted, h2p, "GSW1345", c("bond", controls, adversity))
save_fig(plot_curves(g3, "Predicted probability of adult onset",
  "Figure 3. Strain and peer deviance: adult-onset offending or drug use",
  int_caption("H2P", "C", "Weighted logistic Model C (GSW1345, CLUSTER2); Wave I abstainers, onset window Waves III-V.")),
  "figure3_onset_interaction")

# ---------------------------------------------------------------------------
# FIGURE 4: ONSET BY TYPE
# ---------------------------------------------------------------------------
route_lab <- c(offending_only = "Offending only", substance_only = "Illicit drug use only",
               both = "Both")
h2s <- dat[dat$in_h2s, ]
rr <- function(d, rv, label) {
  data.frame(sample = label, route = unname(route_lab),
             pct = sapply(names(route_lab), function(r) 100 * sum(d$GSW1345[d[[rv]] %in% r]) / sum(d$GSW1345)))
}
r4 <- rbind(rr(h2p, "route_all", "Primary: all adult onset (Waves III-V)"),
            rr(h2s, "route_late", "Secondary: late onset only (Waves IV-V)"))
r4$route <- factor(r4$route, levels = rev(unname(route_lab)))
p4 <- ggplot(r4, aes(pct, route)) + geom_col(fill = "grey45", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.15, size = 3.2, family = "serif") +
  facet_wrap(~ sample, ncol = 1) + scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Weighted % of the at-risk sample", y = NULL,
       title = "Figure 4. Kinds of adult onset among abstainers",
       caption = "Illicit drug use includes marijuana; alcohol excluded. Arrest not used.")
save_fig(p4, "figure4_onset_types", w = 7.5, h = 4.5)

# ---------------------------------------------------------------------------
# FIGURE 5: WEIGHTED VS. UNWEIGHTED
# ---------------------------------------------------------------------------
focal <- c(strain = "Strain", peer = "Peer deviance", bond = "Bond deficit",
           "strain:peer" = "Strain x Peer", phys_abuse = "Physical abuse",
           sex_abuse = "Sexual abuse", adhd = "ADHD symptoms")
spec_lab <- c(H1 = "H1: Wave I delinquency (B)", H2P = "Primary: all adult onset (C)",
              H2S = "Secondary: late onset (C)")
f5 <- ct[ct$family == "logit" & ct$term %in% names(focal) &
         ((ct$spec == "H1" & ct$model == "B") | (ct$spec %in% c("H2P", "H2S") & ct$model == "C")), ]
f5$term_lab <- factor(focal[f5$term], levels = rev(unname(focal)))
f5$spec_lab <- factor(spec_lab[f5$spec], levels = unname(spec_lab))
p5 <- ggplot(f5, aes(OR, term_lab, shape = weighting)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = OR_lo, xmax = OR_hi), position = position_dodge(width = 0.5), size = 0.3) +
  facet_wrap(~ spec_lab, ncol = 3) + scale_x_log10() +
  scale_shape_manual(values = c(weighted = 16, unweighted = 1)) +
  labs(x = "Odds ratio (log scale), 95% CI", y = NULL, shape = NULL,
       title = "Figure 5. Focal and early-adversity effects, weighted (primary) vs. unweighted")
save_fig(p5, "figure5_weighted_vs_unweighted", w = 9.5, h = 4.5)

# ---------------------------------------------------------------------------
# FIGURE 6: OFFENDING ONSET VS. ILLICIT DRUG USE ONSET
# ---------------------------------------------------------------------------
f6 <- ct[ct$family == "logit" & ct$weighting == "weighted" & ct$model == "C" &
         ct$spec %in% c("H2P_OFF", "H2P_SUB", "H2S_OFF", "H2S_SUB") & ct$term %in% names(focal), ]
f6$term_lab <- factor(focal[f6$term], levels = rev(unname(focal)))
f6$design <- ifelse(grepl("^H2P", f6$spec), "Primary: all adult onset (W3-W5)", "Secondary: late onset (W4-W5)")
f6$outcome <- ifelse(grepl("OFF", f6$spec), "Offending onset", "Illicit drug use onset")
f6$design <- paste0(f6$design, "\n(offending ", ave(f6$n_events, f6$design, FUN = function(x) x[1]), " / drug ",
                    ave(f6$n_events, f6$design, FUN = function(x) x[length(x)]), " events)")
p6 <- ggplot(f6, aes(OR, term_lab, shape = outcome)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = OR_lo, xmax = OR_hi), position = position_dodge(width = 0.5), size = 0.3) +
  facet_wrap(~ design) + scale_x_log10() + scale_shape_manual(values = c(17, 16)) +
  labs(x = "Odds ratio (log scale), 95% CI", y = NULL, shape = NULL,
       title = "Figure 6. Predictors of offending onset vs. illicit drug use onset",
       caption = "Weighted Model C (GSW1345). Each outcome estimated on the full risk set.")
save_fig(p6, "figure6_offending_vs_drug_onset", w = 9, h = 4.8)

cat("Figures saved to output/figures/\n")
