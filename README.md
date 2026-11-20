# Traditional and herbal medicine use for children by caregivers in Nigeria: systematic review and meta-analysis

Data and code for the systematic review registered in PROSPERO (CRD420261552693): "Prevalence and determinants of traditional and herbal medicine use for children by caregivers in Nigeria: a systematic review and meta-analysis".

Authors: Olaniyi Felix Sanni and Abike Elizabeth Sanni, Research Department, EasyGlob Health Initiatives, Ota, Ogun State, Nigeria.

## Contents

- `data/meta_analysis_data.csv`: the 134 studies in the pooled prevalence analysis (events, denominators, condition, zone, setting, period, type of traditional medicine, risk of bias, recall period, child age and denominator unit groups)
- `data/included_studies_characteristics.csv`: characteristics of the 159 included reports of primary studies and secondary analyses
- `data/primary_estimates.csv`: the primary estimate chosen for each study and whether it was pooled
- `data/all_tm_measures.csv`: every measure of traditional medicine use extracted (423 measures)
- `data/associated_factors.csv`: associations with traditional medicine use as extracted
- `data/conditions_remedies_safety.csv`: conditions treated, remedy types, concurrent use, disclosure and adverse events
- `data/risk_of_bias_jbi.csv`: Joanna Briggs Institute prevalence checklist ratings (Q1 to Q9) for each report
- `data/national_survey_reports.csv`: assessment and extracted figures for the 11 NDHS and MICS reports
- `data/fulltext_screening_decisions.csv`: full-text decisions with reasons for all 329 reports assessed
- `data/S3_excluded_reports.csv`: reports excluded at full text, with reasons
- `data/prisma_counts.json`: record counts for the PRISMA 2020 flow diagram
- `code/run_meta_analysis.R`: all pooled, subgroup, sensitivity, leave-one-out and meta-regression analyses, plus the exploratory subgroups and joint model
- `results/`: output of the script, and the synthesis-without-meta-analysis tables (factors, conditions, remedies, disclosure and adverse events)
- `figures/`: forest and funnel plots produced on MetaAnalysisOnline.com, the PRISMA flow diagram and the risk of bias figure
- `search/search_strategies.md`: full search strategies, dates and record counts

## How to run

From the repository root, with R 4.3.3, meta 9.0.0 and metafor 4.4-0:

```
Rscript code/run_meta_analysis.R
```

Results are written to `results/` and the R check figures to `results/r_figures/`.

## Methods summary

Proportions were pooled with a random-effects generalized linear mixed model (logit link, maximum likelihood), with 95% confidence and prediction intervals. Inverse-variance pooling of logit proportions with REML and the Knapp-Hartung adjustment (MetaAnalysisOnline.com) was the comparison model and produced the forest plots. Subgroup differences were tested with mixed-effects logistic meta-regression assuming a common between-study variance. National survey reports were tabulated as a separate stratum. Associated factors were synthesized by direction of effect with sign tests.

## Licence

Data: CC BY 4.0. Code: MIT.
