# Review 4: Traditional and herbal medicine use for children by caregivers in Nigeria
# Meta-analysis of proportions. Primary model: GLMM (random-intercept logistic), logit scale.
# Comparison model: inverse variance on logit, REML tau2, Hartung-Knapp CI (IQWiG6 ad hoc correction).
suppressPackageStartupMessages({library(meta); library(metafor)})
# Run from the repository root: Rscript code/run_meta_analysis.R
dir.create("results/r_figures", recursive = TRUE, showWarnings = FALSE)
d <- read.csv("data/meta_analysis_data.csv", stringsAsFactors = FALSE)

# GoogleScholar_0031 reports only 3.4% (95% CI 2.73-4.22): effective n and events from the CI on the logit scale
i31 <- which(d$rec_id == "GoogleScholar_0031")
if (length(i31) == 1) {
  p <- 0.034; se <- (qlogis(0.0422) - qlogis(0.0273)) / 3.92
  neff <- round(1 / (se^2 * p * (1 - p))); d$n[i31] <- neff; d$event[i31] <- round(p * neff)
}
d$event <- as.numeric(d$event); d$n <- as.numeric(d$n)
d$year <- suppressWarnings(as.numeric(sub(".*?(\\d{4}).*", "\\1", d$study)))
cat("Studies:", nrow(d), " Participants:", sum(d$n), "\n")

fit <- function(dat, sub = NULL, ...) {
  metaprop(event, n, studlab = study, data = dat, sm = "PLOGIT", method = "GLMM",
           method.tau = "ML", common = FALSE, random = TRUE, prediction = TRUE,
           subgroup = sub, ...)
}
row_of <- function(m, label) {
  data.frame(analysis = label, k = m$k, n = sum(m$n),
             pooled = round(100 * plogis(m$TE.random), 1),
             lcl = round(100 * plogis(m$lower.random), 1),
             ucl = round(100 * plogis(m$upper.random), 1),
             pi_l = round(100 * plogis(m$lower.predict), 1),
             pi_u = round(100 * plogis(m$upper.predict), 1),
             I2 = round(100 * m$I2, 1), tau2 = round(m$tau2, 3))
}

# 1. Primary GLMM
m0 <- fit(d)
res <- row_of(m0, "Primary: GLMM, all poolable studies")

# 2. Comparison: IV REML + Hartung-Knapp
m_iv <- metaprop(event, n, studlab = study, data = d, sm = "PLOGIT", method = "Inverse",
                 method.tau = "REML", method.random.ci = "HK", adhoc.hakn.ci = "IQWiG6",
                 common = FALSE, prediction = TRUE, incr = 0.5)
res <- rbind(res, row_of(m_iv, "Comparison: inverse variance, REML, Hartung-Knapp"))
m_dl <- update(m_iv, method.tau = "DL", method.random.ci = "classic")
res <- rbind(res, row_of(m_dl, "Comparison: inverse variance, DerSimonian-Laird"))

# 3. Subgroups (GLMM), test for subgroup differences
subs <- c(zone = "Geopolitical zone", condition = "Condition", setting = "Setting",
          period = "Publication period", tm_type = "Type of traditional medicine",
          rob = "Risk of bias (JBI)", respondent = "Respondent")
sg <- list(); sgtest <- data.frame()
for (v in names(subs)) {
  dd <- d[!is.na(d[[v]]) & d[[v]] != "", ]
  m <- fit(dd, sub = dd[[v]])
  sg[[v]] <- m
  lev <- m$subgroup.levels
  for (j in seq_along(lev)) {
    sgtest <- rbind(sgtest, data.frame(variable = subs[v], level = lev[j], k = m$k.w[j],
      pooled = round(100 * plogis(m$TE.random.w[j]), 1),
      lcl = round(100 * plogis(m$lower.random.w[j]), 1),
      ucl = round(100 * plogis(m$upper.random.w[j]), 1),
      I2 = round(100 * m$I2.w[j], 1),
      p_between = signif(m$pval.Q.b.random, 3)))
  }
}
write.csv(sgtest, "results/subgroups.csv", row.names = FALSE)

# 4. Sensitivity analyses
sens <- list(
  "Excluding high risk of bias" = d$rob != "high",
  "Low risk of bias only" = d$rob == "low",
  "Excluding flagged studies (implausible data, age 18, preventive use)" = !grepl("implausible|age18|preventive", d$flags),
  "Excluding denominators under 50" = d$n >= 50,
  "Community and national surveys only" = d$setting %in% c("Community/household", "National survey"),
  "Excluding GoogleScholar_0031 (CI-derived counts)" = d$rec_id != "GoogleScholar_0031"
)
for (s in names(sens)) res <- rbind(res, row_of(fit(d[sens[[s]], ]), s))
write.csv(res, "results/pooled_and_sensitivity.csv", row.names = FALSE)

# 5. Leave-one-out and influence (IV REML model, logit)
es <- escalc(measure = "PLO", xi = event, ni = n, data = d, add = 0.5, to = "only0")
r_reml <- rma(yi, vi, data = es, method = "REML", test = "knha")
loo <- leave1out(r_reml)
loo_df <- data.frame(omitted = d$study, pooled = round(100 * plogis(loo$estimate), 1),
                     lcl = round(100 * plogis(loo$ci.lb), 1), ucl = round(100 * plogis(loo$ci.ub), 1),
                     I2 = round(loo$I2, 1))
write.csv(loo_df, "results/leave_one_out.csv", row.names = FALSE)
cat("Leave-one-out range:", range(loo_df$pooled), "\n")

# 6. Meta-regression (GLMM, logit) on publication year and zone
mr_year <- rma.glmm(measure = "PLO", xi = event, ni = n, mods = ~ year, data = d[!is.na(d$year), ])
mr_zone <- rma.glmm(measure = "PLO", xi = event, ni = n, mods = ~ factor(zone), data = d)
mr_multi <- rma.glmm(measure = "PLO", xi = event, ni = n,
                     mods = ~ year + factor(setting) + factor(rob), data = d[!is.na(d$year), ])
sink("results/meta_regression.txt")
cat("Meta-regression on publication year\n"); print(mr_year)
cat("\nMeta-regression on zone\n"); print(mr_zone)
cat("\nMultivariable: year + setting + risk of bias\n"); print(mr_multi)
sink()

# 7. Small-study effects: funnel plot (descriptive only for prevalence; Egger shown for transparency, not inference)
png("results/r_figures/fig_funnel_primary.png", width = 7, height = 6, units = "in", res = 600)
funnel(r_reml, xlab = "Logit proportion", main = "Funnel plot, logit proportion (descriptive)")
dev.off()
eg <- regtest(r_reml)
sink("results/small_study.txt"); print(eg); sink()

# 8. Forest plots
h <- function(k) max(8, 0.17 * k + 3)
png("results/r_figures/fig_forest_primary_by_condition.png", width = 11, height = h(nrow(d) + 30), units = "in", res = 300)
forest(sg$condition, sortvar = TE, leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Users", "Total"), rightcols = c("effect", "ci"),
       rightlabs = c("Proportion", "95% CI"), pscale = 100, digits = 1,
       print.subgroup.name = FALSE, col.diamond = "navy", fontsize = 8, spacing = 0.9,
       xlab = "Proportion using traditional or herbal medicine for the child (%)",
       prediction = TRUE, test.subgroup = TRUE)
dev.off()
pdf("results/r_figures/fig_forest_primary_by_condition.pdf", width = 11, height = h(nrow(d) + 30))
forest(sg$condition, sortvar = TE, leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Users", "Total"), pscale = 100, digits = 1,
       print.subgroup.name = FALSE, fontsize = 8, spacing = 0.9, prediction = TRUE, test.subgroup = TRUE)
dev.off()
for (v in c("zone", "setting", "period", "tm_type", "rob")) {
  png(sprintf("results/r_figures/fig_forest_subgroup_%s.png", v), width = 9, height = 4 + 0.5 * length(sg[[v]]$subgroup.levels), units = "in", res = 600)
  forest(sg[[v]], study.results = FALSE, pscale = 100, digits = 1, leftcols = "studlab", leftlabs = subs[v],
         print.subgroup.name = FALSE, test.subgroup = TRUE, prediction = FALSE, fontsize = 10,
         xlab = "Pooled proportion (%)")
  dev.off()
}
saveRDS(list(m0 = m0, m_iv = m_iv, sg = sg, res = res), "results/models.rds")
print(res)
print(sgtest)
sink("results/sessionInfo.txt"); print(sessionInfo()); sink()

# 9. Subgroup differences with a common between-study variance (GLMM meta-regression, omnibus QM test)
qm <- data.frame()
for (v in names(subs)) {
  dd <- d[!is.na(d[[v]]) & d[[v]] != "", ]
  r <- try(rma.glmm(measure = "PLO", xi = event, ni = n, mods = as.formula(paste0("~ factor(", v, ")")), data = dd), silent = TRUE)
  if (!inherits(r, "try-error")) qm <- rbind(qm, data.frame(variable = subs[v], QM = round(r$QM, 2), df = r$QMdf[1], p = signif(r$QMp, 3), tau2 = round(r$tau2, 3)))
}
# zone without national/multi-state studies
dz <- d[!d$zone %in% c("National", "Multiple"), ]
rz <- rma.glmm(measure = "PLO", xi = event, ni = n, mods = ~ factor(zone), data = dz)
qm <- rbind(qm, data.frame(variable = "Geopolitical zone (single-zone studies only)", QM = round(rz$QM, 2), df = rz$QMdf[1], p = signif(rz$QMp, 3), tau2 = round(rz$tau2, 3)))
write.csv(qm, "results/subgroup_tests_common_tau.csv", row.names = FALSE)
print(qm)


# 10. Exploratory subgroups: recall period, child age, denominator unit (common tau2 tests)
out<-data.frame(); qm<-data.frame()
for(v in c("age_group","recall_group","unit_group")){
  m<-metaprop(event,n,data=d,sm="PLOGIT",method="GLMM",method.tau="ML",subgroup=d[[v]],common=FALSE)
  for(j in seq_along(m$subgroup.levels)) out<-rbind(out,data.frame(variable=v,level=m$subgroup.levels[j],k=m$k.w[j],n=sum(d$n[d[[v]]==m$subgroup.levels[j]]),pooled=round(100*plogis(m$TE.random.w[j]),1),lcl=round(100*plogis(m$lower.random.w[j]),1),ucl=round(100*plogis(m$upper.random.w[j]),1)))
  r<-rma.glmm(measure="PLO",xi=event,ni=n,mods=as.formula(paste0("~factor(",v,")")),data=d)
  qm<-rbind(qm,data.frame(variable=v,QM=round(r$QM,2),df=r$QMdf[1],p=signif(r$QMp,3),tau2=round(r$tau2,3)))
}
print(out);print(qm)
write.csv(out,"results/subgroups_extra.csv",row.names=FALSE);write.csv(qm,"results/subgroup_tests_extra.csv",row.names=FALSE)


# 11. Sensitivity: excluding secondary analyses of national surveys
dd<-d[d$setting!="National survey",]
m<-metaprop(event,n,data=dd,sm="PLOGIT",method="GLMM",method.tau="ML",prediction=TRUE,common=FALSE)
cat(nrow(dd),sum(dd$n),round(100*plogis(c(m$TE.random,m$lower.random,m$upper.random,m$lower.predict,m$upper.predict)),1),round(m$I2*100,1),round(m$tau2,3),"\n")
# recall x condition quick look
print(table(d$recall_group,d$condition))


# 12. Joint model (recall period + condition), case-control sensitivity and sign-test checks
d2<-d; d2$recall_group<-relevel(factor(d$recall_group),ref="Specific illness episode")
r<-rma.glmm(measure="PLO",xi=event,ni=n,mods=~recall_group+factor(condition),data=d2,control=list(optimizer="bobyqa"))
print(r)
b<-coef(r);V<-vcov(r);f<-function(ix){w<-t(b[ix])%*%solve(V[ix,ix])%*%b[ix];c(w,length(ix),pchisq(w,length(ix),lower.tail=FALSE))}
cat("tau2",r$tau2,"QM",r$QM,r$QMp,"\ncondition Wald",f(4:10),"\nrecall Wald",f(2:3),"\nOR usual vs episode",exp(c(b[3],b[3]-1.96*sqrt(V[3,3]),b[3]+1.96*sqrt(V[3,3]))),"\n")
m<-metaprop(event,n,data=d[d$rec_id!="AJOL_0031",],sm="PLOGIT",method="GLMM",method.tau="ML",prediction=TRUE,common=FALSE)
cat("excl Ibadin:",sum(d$rec_id!="AJOL_0031"),sum(d$n[d$rec_id!="AJOL_0031"]),round(100*plogis(c(m$TE.random,m$lower.random,m$upper.random,m$lower.predict,m$upper.predict)),1),round(100*m$I2,1),round(m$tau2,2),"\n")
cat(binom.test(15,16)$p.value, binom.test(8,9)$p.value, binom.test(15,20)$p.value, binom.test(8,11)$p.value,"\n")
